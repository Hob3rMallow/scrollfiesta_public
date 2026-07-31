/* Mesh quality analyzer: distribution of triangle areas, edge lengths,
 * aspect ratios, connected components, boundary loop sizes. */
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { float z, y, x; } V3;

static V3 *V; static int nv = 0; static int vcap = 1<<16;
static int *F; static int nf = 0; static int fcap = 1<<16;

static int *uf_parent;
static int uf_find(int i) { while (uf_parent[i] != i) { uf_parent[i] = uf_parent[uf_parent[i]]; i = uf_parent[i]; } return i; }
static void uf_union(int a, int b) { a=uf_find(a); b=uf_find(b); if (a!=b) uf_parent[a]=b; }

int main(int argc, char **argv) {
    if (argc<2) { fprintf(stderr,"usage: %s mesh.obj\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1],"r"); if (!f) { perror(argv[1]); return 1; }
    V = malloc(vcap*sizeof(V3));
    F = malloc(fcap*3*sizeof(int));
    char line[256];
    while (fgets(line,sizeof(line),f)) {
        if (line[0]=='v' && line[1]==' ') {
            V3 v;
            if (sscanf(line, "v %f %f %f", &v.z, &v.y, &v.x) == 3) {
                if (nv>=vcap) { vcap*=2; V=realloc(V,vcap*sizeof(V3)); }
                V[nv++] = v;
            }
        } else if (line[0]=='f' && line[1]==' ') {
            int a,b,c;
            if (sscanf(line,"f %d %d %d",&a,&b,&c)==3) {
                if (nf>=fcap) { fcap*=2; F=realloc(F,fcap*3*sizeof(int)); }
                F[nf*3+0]=a-1; F[nf*3+1]=b-1; F[nf*3+2]=c-1; nf++;
            }
        }
    }
    fclose(f);
    fprintf(stderr,"Loaded %d verts, %d faces\n", nv, nf);

    /* Triangle area + edge length stats */
    double area_min=1e30, area_max=0, area_sum=0;
    double edge_min=1e30, edge_max=0, edge_sum=0; int edge_n=0;
    double ar_max=0; double ar_sum=0;
    int n_degenerate = 0;
    int area_hist[20] = {0};  /* log-binned: 1e-6 .. 1e6 */
    for (int i=0; i<nf; i++) {
        int a=F[i*3+0], b=F[i*3+1], c=F[i*3+2];
        if (a<0||a>=nv||b<0||b>=nv||c<0||c>=nv) { n_degenerate++; continue; }
        double az=V[a].z, ay=V[a].y, ax=V[a].x;
        double bz=V[b].z, by=V[b].y, bx=V[b].x;
        double cz=V[c].z, cy=V[c].y, cx=V[c].x;
        double e1z=bz-az, e1y=by-ay, e1x=bx-ax;
        double e2z=cz-az, e2y=cy-ay, e2x=cx-ax;
        double e3z=cz-bz, e3y=cy-by, e3x=cx-bx;
        double l1=sqrt(e1z*e1z+e1y*e1y+e1x*e1x);
        double l2=sqrt(e2z*e2z+e2y*e2y+e2x*e2x);
        double l3=sqrt(e3z*e3z+e3y*e3y+e3x*e3x);
        double nz=e1y*e2x-e1x*e2y, ny=e1x*e2z-e1z*e2x, nx=e1z*e2y-e1y*e2z;
        double area = 0.5 * sqrt(nz*nz+ny*ny+nx*nx);
        if (area<1e-12) { n_degenerate++; continue; }
        area_sum += area;
        if (area<area_min) area_min=area;
        if (area>area_max) area_max=area;
        edge_sum += l1+l2+l3; edge_n+=3;
        if (l1<edge_min) edge_min=l1; if (l2<edge_min) edge_min=l2; if (l3<edge_min) edge_min=l3;
        if (l1>edge_max) edge_max=l1; if (l2>edge_max) edge_max=l2; if (l3>edge_max) edge_max=l3;
        double lmax=l1; if (l2>lmax) lmax=l2; if (l3>lmax) lmax=l3;
        double lmin=l1; if (l2<lmin) lmin=l2; if (l3<lmin) lmin=l3;
        double ar = lmax/lmin;
        ar_sum+=ar; if (ar>ar_max) ar_max=ar;
        int bin = (int)(log10(area*1e6)+0.5); if (bin<0) bin=0; if (bin>=20) bin=19;
        area_hist[bin]++;
    }

    fprintf(stderr,"\n=== Triangle areas ===\n");
    fprintf(stderr,"  min=%.4g  max=%.4g  mean=%.4g  degenerate=%d\n",
            area_min, area_max, area_sum/(nf-n_degenerate), n_degenerate);
    fprintf(stderr,"  log-binned (area*1e6) histogram:\n");
    for (int i=0; i<20; i++)
        if (area_hist[i]>0) fprintf(stderr,"    1e%d..1e%d: %d\n", i-6, i-5, area_hist[i]);

    fprintf(stderr,"\n=== Edge lengths ===\n");
    fprintf(stderr,"  min=%.4f  max=%.4f  mean=%.4f vox\n",
            edge_min, edge_max, edge_sum/edge_n);

    fprintf(stderr,"\n=== Aspect ratios (long-edge / short-edge) ===\n");
    fprintf(stderr,"  max=%.2f  mean=%.2f\n", ar_max, ar_sum/(nf-n_degenerate));

    /* Connected components by face graph */
    uf_parent = malloc(nv*sizeof(int));
    for (int i=0; i<nv; i++) uf_parent[i]=i;
    for (int i=0; i<nf; i++) {
        int a=F[i*3+0], b=F[i*3+1], c=F[i*3+2];
        if (a<0||a>=nv||b<0||b>=nv||c<0||c>=nv) continue;
        uf_union(a,b); uf_union(b,c);
    }
    /* Count CC sizes */
    int *cc_size = calloc(nv, sizeof(int));
    int used_count = 0;
    int *vert_used = calloc(nv, sizeof(int));
    for (int i=0; i<nf; i++) {
        int a=F[i*3+0], b=F[i*3+1], c=F[i*3+2];
        if (a<0||a>=nv||b<0||b>=nv||c<0||c>=nv) continue;
        vert_used[a] = vert_used[b] = vert_used[c] = 1;
    }
    for (int i=0; i<nv; i++) if (vert_used[i]) { cc_size[uf_find(i)]++; used_count++; }

    int ncc = 0;
    int large_cc = 0, small_cc = 0;
    int sizes[100] = {0};
    for (int i=0; i<nv; i++) {
        if (cc_size[i] > 0) {
            ncc++;
            if (cc_size[i] > 100) large_cc++; else small_cc++;
            int bin = cc_size[i]; if (bin >= 100) bin = 99;
            sizes[bin]++;
        }
    }
    fprintf(stderr,"\n=== Connected components (on the face graph) ===\n");
    fprintf(stderr,"  total CCs: %d  (large >100v: %d, small ≤100v: %d)\n",
            ncc, large_cc, small_cc);
    fprintf(stderr,"  used verts: %d / %d (%.1f%%)\n",
            used_count, nv, 100.0*used_count/nv);
    fprintf(stderr,"  small-CC size histogram (≤30 verts):\n");
    for (int i=3; i<=30 && i<100; i++)
        if (sizes[i]) fprintf(stderr,"    size=%d: %d CCs\n", i, sizes[i]);
    /* Histogram of larger sizes via binning */
    int large_bins[10] = {0};
    for (int i=0; i<nv; i++) {
        if (cc_size[i] > 100) {
            int b = (int)log10(cc_size[i]); if (b>=10) b=9;
            large_bins[b]++;
        }
    }
    fprintf(stderr,"  large-CC size histogram (decimal log of vert count):\n");
    for (int i=2; i<10; i++)
        if (large_bins[i]) fprintf(stderr,"    1e%d-1e%d: %d CCs\n", i, i+1, large_bins[i]);

    free(uf_parent); free(cc_size); free(vert_used);
    free(V); free(F);
    return 0;
}
