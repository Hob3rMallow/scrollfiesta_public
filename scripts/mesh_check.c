/* Manifoldness check on triangle mesh OBJ. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define EDGE_HASH_SIZE (1u<<22)

typedef struct EdgeKey { int va, vb, count, next; } EdgeKey;

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s mesh.obj\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }
    int cap=1<<16, nv=0, nf=0;
    int *F = (int *)malloc(cap*3*sizeof(int));
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0]=='v' && line[1]==' ') nv++;
        else if (line[0]=='f' && line[1]==' ') {
            int a,b,c;
            if (sscanf(line, "f %d %d %d", &a, &b, &c) == 3) {
                if (nf>=cap) { cap*=2; F=realloc(F, cap*3*sizeof(int)); }
                F[nf*3+0]=a-1; F[nf*3+1]=b-1; F[nf*3+2]=c-1;
                nf++;
            }
        }
    }
    fclose(f);
    int *bk = (int *)malloc(EDGE_HASH_SIZE*sizeof(int));
    for (size_t i=0; i<EDGE_HASH_SIZE; i++) bk[i]=-1;
    int key_cap = nf*3, nk = 0;
    EdgeKey *K = (EdgeKey *)calloc(key_cap, sizeof(EdgeKey));
    for (int i=0; i<nf; i++) {
        int tri[3] = {F[i*3], F[i*3+1], F[i*3+2]};
        for (int e=0; e<3; e++) {
            int va=tri[e], vb=tri[(e+1)%3];
            if (va>vb) { int t=va; va=vb; vb=t; }
            uint64_t h = (uint64_t)(uint32_t)va * 0x9E3779B97F4A7C15ULL + vb;
            h ^= h>>33; h *= 0xff51afd7ed558ccdULL; h ^= h>>33;
            size_t b = h & (EDGE_HASH_SIZE-1);
            int idx = bk[b];
            while (idx>=0) {
                if (K[idx].va==va && K[idx].vb==vb) break;
                idx = K[idx].next;
            }
            if (idx<0) {
                if (nk >= key_cap) { key_cap *= 2; K = realloc(K, key_cap*sizeof(EdgeKey)); }
                idx = nk++;
                K[idx].va=va; K[idx].vb=vb; K[idx].count=0;
                K[idx].next = bk[b]; bk[b] = idx;
            }
            K[idx].count++;
        }
    }
    int boundary=0, manifold=0, nonmanifold=0;
    for (int i=0; i<nk; i++) {
        if (K[i].count==1) boundary++;
        else if (K[i].count==2) manifold++;
        else nonmanifold++;
    }
    int *used = (int *)calloc(nv, sizeof(int));
    for (int i=0; i<nf; i++) {
        used[F[i*3]]++; used[F[i*3+1]]++; used[F[i*3+2]]++;
    }
    int isolated=0;
    for (int i=0; i<nv; i++) if (!used[i]) isolated++;
    printf("Mesh: %d verts, %d faces, %d unique edges\n", nv, nf, nk);
    printf("  boundary edges (1 face):       %d\n", boundary);
    printf("  manifold interior (2 faces):   %d\n", manifold);
    printf("  NON-manifold (3+ faces):       %d\n", nonmanifold);
    printf("  isolated verts:                %d\n", isolated);
    printf("  V - E + F = %d\n", nv - nk + nf);
    free(F); free(bk); free(K); free(used);
    return 0;
}
