/* Filter an OBJ point cloud (v + vn) to keep only verts inside a world bbox. */
#include <stdio.h>
#include <stdlib.h>
#define MAX_VERTS 1000000
typedef struct { float z,y,x; } V3;
int main(int argc, char **argv)
{
    if (argc < 9) { fprintf(stderr, "usage: %s in.obj out.obj zlo zhi ylo yhi xlo xhi\n", argv[0]); return 1; }
    float zlo=atof(argv[3]), zhi=atof(argv[4]), ylo=atof(argv[5]), yhi=atof(argv[6]), xlo=atof(argv[7]), xhi=atof(argv[8]);
    FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    V3 *V = malloc(MAX_VERTS*sizeof(V3)), *N = malloc(MAX_VERTS*sizeof(V3));
    int nv=0, nvn=0; char line[256];
    while (fgets(line, sizeof(line), f)) {
        V3 v;
        if (line[0]=='v' && line[1]==' ') {
            if (sscanf(line, "v %f %f %f", &v.z, &v.y, &v.x)==3 && nv<MAX_VERTS) V[nv++]=v;
        } else if (line[0]=='v' && line[1]=='n' && line[2]==' ') {
            if (sscanf(line, "vn %f %f %f", &v.z, &v.y, &v.x)==3 && nvn<MAX_VERTS) N[nvn++]=v;
        }
    }
    fclose(f);
    FILE *o = fopen(argv[2], "w"); if (!o) { perror(argv[2]); return 1; }
    int *keep = calloc(nv, sizeof(int)); int kept=0;
    for (int i=0; i<nv; i++)
        if (V[i].z>=zlo && V[i].z<=zhi && V[i].y>=ylo && V[i].y<=yhi && V[i].x>=xlo && V[i].x<=xhi) { keep[i]=1; kept++; }
    for (int i=0; i<nv; i++) if (keep[i]) fprintf(o, "v %.6f %.6f %.6f\n", V[i].z, V[i].y, V[i].x);
    for (int i=0; i<nvn; i++) if (keep[i]) fprintf(o, "vn %.6f %.6f %.6f\n", N[i].z, N[i].y, N[i].x);
    fclose(o);
    fprintf(stderr, "kept %d/%d\n", kept, nv);
    return 0;
}
