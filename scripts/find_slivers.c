/* Find the worst sliver triangles in an OBJ mesh and print details. */
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { float z, y, x; } V3;

int main(int argc, char **argv) {
    if (argc<2) { fprintf(stderr,"usage: %s mesh.obj [n_worst=10]\n", argv[0]); return 1; }
    int n_worst = (argc>=3)?atoi(argv[2]):10;
    FILE *f=fopen(argv[1],"r"); if (!f) { perror(argv[1]); return 1; }
    int vcap=1<<16, nv=0, fcap=1<<16, nf=0;
    V3 *V=malloc(vcap*sizeof(V3));
    int *F=malloc(fcap*3*sizeof(int));
    char line[256];
    while (fgets(line,sizeof(line),f)) {
        if (line[0]=='v' && line[1]==' ') {
            V3 v;
            if (sscanf(line, "v %f %f %f", &v.z, &v.y, &v.x) == 3) {
                if (nv>=vcap){vcap*=2;V=realloc(V,vcap*sizeof(V3));}
                V[nv++]=v;
            }
        } else if (line[0]=='f' && line[1]==' ') {
            int a,b,c;
            if (sscanf(line,"f %d %d %d",&a,&b,&c)==3) {
                if (nf>=fcap){fcap*=2;F=realloc(F,fcap*3*sizeof(int));}
                F[nf*3+0]=a-1; F[nf*3+1]=b-1; F[nf*3+2]=c-1; nf++;
            }
        }
    }
    fclose(f);

    /* For each triangle compute (aspect_ratio, min_edge, area) */
    typedef struct { int fi; double ar; double min_edge; double area; } Q;
    Q *qs = malloc(nf*sizeof(Q));
    for (int i=0; i<nf; i++) {
        int a=F[i*3], b=F[i*3+1], c=F[i*3+2];
        double az=V[a].z, ay=V[a].y, ax=V[a].x;
        double bz=V[b].z, by=V[b].y, bx=V[b].x;
        double cz=V[c].z, cy=V[c].y, cx=V[c].x;
        double e1=sqrt((bz-az)*(bz-az)+(by-ay)*(by-ay)+(bx-ax)*(bx-ax));
        double e2=sqrt((cz-az)*(cz-az)+(cy-ay)*(cy-ay)+(cx-ax)*(cx-ax));
        double e3=sqrt((cz-bz)*(cz-bz)+(cy-by)*(cy-by)+(cx-bx)*(cx-bx));
        double lmin=e1; if(e2<lmin)lmin=e2; if(e3<lmin)lmin=e3;
        double lmax=e1; if(e2>lmax)lmax=e2; if(e3>lmax)lmax=e3;
        double nz=(by-ay)*(cx-ax)-(bx-ax)*(cy-ay);
        double ny=(bx-ax)*(cz-az)-(bz-az)*(cx-ax);
        double nx=(bz-az)*(cy-ay)-(by-ay)*(cz-az);
        double area=0.5*sqrt(nz*nz+ny*ny+nx*nx);
        qs[i].fi=i; qs[i].ar = (lmin>0)?lmax/lmin:1e30;
        qs[i].min_edge = lmin;
        qs[i].area = area;
    }
    /* sort by aspect ratio descending */
    /* simple n*log(n) qsort would do; use a heap-extract for top-N */
    for (int i=0; i<n_worst; i++) {
        int best=i;
        for (int j=i+1; j<nf; j++) if (qs[j].ar > qs[best].ar) best=j;
        Q tmp=qs[i]; qs[i]=qs[best]; qs[best]=tmp;
    }
    printf("Top %d worst-aspect triangles:\n", n_worst);
    printf("  rank  face_idx  aspect   min_edge   area      verts (z,y,x)\n");
    for (int i=0; i<n_worst; i++) {
        Q q=qs[i];
        int a=F[q.fi*3], b=F[q.fi*3+1], c=F[q.fi*3+2];
        printf("  %2d    %7d   %7.2f  %.5f   %.5f   v%d=(%.3f,%.3f,%.3f) v%d=(%.3f,%.3f,%.3f) v%d=(%.3f,%.3f,%.3f)\n",
            i+1, q.fi, q.ar, q.min_edge, q.area,
            a, V[a].z, V[a].y, V[a].x,
            b, V[b].z, V[b].y, V[b].x,
            c, V[c].z, V[c].y, V[c].x);
    }
    free(qs); free(V); free(F);
    return 0;
}
