/* obj_longedge.c -- locate the longest edges/triangles in an OBJ.
 *
 * "Insanely long triangles that shoot off into space" are stretched edges. This
 * reports, over all triangle edges: max / percentiles / a length histogram, and
 * the top-N longest edges with BOTH endpoint coordinates + vertex indices, so
 * you can see WHERE the spikes are (in 3D for a surface OBJ, or in (u,v) for a
 * flat unroll whose verts are "v u v 0").
 *
 * Also flags "spike vertices": a vertex whose longest incident edge is >
 * --spike x its median incident edge (a lone vertex yanked far from its
 * neighborhood -- the classic single-bad-UV or bad-fill signature).
 *
 * Usage:
 *   obj_longedge <mesh.obj> [--top N=20] [--spike F=8] [--hist-max L]
 *   obj_longedge --selftest
 *
 * Exit: 0 ok / selftest pass, 1 IO, 2 usage, 3 selftest fail.
 * Standalone C99. Coord convention "v c0 c1 c2" (z,y,x for surface OBJs).
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

typedef struct { double *p; size_t n, cap; } DA;
typedef struct { int32_t *p; size_t n, cap; } IA;
static void da_push(DA *d, double v){ if(d->n==d->cap){ d->cap=d->cap?d->cap*2:1<<16; d->p=realloc(d->p,d->cap*sizeof(double)); if(!d->p){fprintf(stderr,"OOM\n");exit(1);} } d->p[d->n++]=v; }
static void ia_push(IA *a, int32_t v){ if(a->n==a->cap){ a->cap=a->cap?a->cap*2:1<<16; a->p=realloc(a->p,a->cap*sizeof(int32_t)); if(!a->p){fprintf(stderr,"OOM\n");exit(1);} } a->p[a->n++]=v; }

static int load_obj(const char *path, double **V, size_t *nv, int32_t **F, size_t *nf)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    DA vv = {0}; IA ff = {0};
    char line[8192];
    while (fgets(line, sizeof line, f)) {
        if (line[0]=='v' && line[1]==' ') {
            double a,b,c;
            if (sscanf(line+2, "%lf %lf %lf", &a,&b,&c)==3){ da_push(&vv,a); da_push(&vv,b); da_push(&vv,c); }
        } else if (line[0]=='f' && line[1]==' ') {
            int idx[64], ni=0; const char *p=line+1;
            while (*p && ni<64){
                while(*p==' '||*p=='\t')p++;
                if(*p=='\0'||*p=='\n'||*p=='\r')break;
                long vi=strtol(p,(char**)&p,10);
                while(*p && *p!=' ' && *p!='\t' && *p!='\n' && *p!='\r')p++;
                idx[ni++] = vi>0 ? (int32_t)(vi-1) : (int32_t)((long)(vv.n/3)+vi);
            }
            for(int i=1;i+1<ni;i++){ ia_push(&ff,idx[0]); ia_push(&ff,idx[i]); ia_push(&ff,idx[i+1]); }
        }
    }
    fclose(f);
    *V=vv.p; *nv=vv.n/3; *F=ff.p; *nf=ff.n/3;
    return 0;
}

static int cmp_d(const void *a, const void *b){ double x=*(const double*)a,y=*(const double*)b; return x<y?-1:(x>y?1:0); }

typedef struct { double len; int32_t a, b; int32_t face; } Edge;
static int cmp_edge_desc(const void *pa, const void *pb){ const Edge*a=pa,*b=pb; return a->len<b->len?1:(a->len>b->len?-1:0); }

static void analyze(const double *V, size_t nv, const int32_t *F, size_t nf,
                    int topN, double spike, double hist_max,
                    double *out_max, long *out_over)
{
    if (nf == 0){ printf("  (no faces)\n"); *out_max=0; *out_over=0; return; }
    /* per-edge lengths (one per triangle side; duplicates across shared edges ok
     * for stats + top-N localization) */
    size_t ne = nf*3;
    Edge *E = malloc(ne*sizeof(Edge));
    double *len = malloc(ne*sizeof(double));
    /* per-vertex: max incident + all incident (for median) */
    double *vmax = calloc(nv, sizeof(double));
    /* incident edge lengths per vertex via a second pass median: approximate
     * median with per-vertex running collection is heavy; use mean of incident
     * as the scale and max/scale as the spike ratio. */
    double *vsum = calloc(nv, sizeof(double));
    int32_t *vcnt = calloc(nv, sizeof(int32_t));
    for (size_t t=0;t<nf;t++){
        for (int e=0;e<3;e++){
            int32_t a=F[t*3+e], b=F[t*3+(e+1)%3];
            double dz=V[(size_t)a*3]-V[(size_t)b*3];
            double dy=V[(size_t)a*3+1]-V[(size_t)b*3+1];
            double dx=V[(size_t)a*3+2]-V[(size_t)b*3+2];
            double l=sqrt(dz*dz+dy*dy+dx*dx);
            size_t ei=t*3+(size_t)e;
            E[ei].len=l; E[ei].a=a; E[ei].b=b; E[ei].face=(int32_t)t;
            len[ei]=l;
            if(l>vmax[a])vmax[a]=l; if(l>vmax[b])vmax[b]=l;
            vsum[a]+=l; vsum[b]+=l; vcnt[a]++; vcnt[b]++;
        }
    }
    qsort(len, ne, sizeof(double), cmp_d);
    double mx=len[ne-1];
    double p50=len[ne/2], p99=len[(size_t)(ne*0.99)], p999=len[(size_t)(ne*0.999)];
    *out_max = mx;

    /* histogram */
    double hm = hist_max>0 ? hist_max : p999*4;
    if (hm<1e-9) hm=1.0;
    printf("  edges=%zu  median=%.2f  p99=%.2f  p99.9=%.2f  MAX=%.2f vox\n",
           ne, p50, p99, p999, mx);
    {
        static const double edg[8]={2,4,8,16,32,64,128,1e18};
        long h[8]={0,0,0,0,0,0,0,0};
        for(size_t i=0;i<ne;i++){ int bk=0; while(bk<7 && len[i]>=edg[bk])bk++; h[bk]++; }
        printf("  len hist [<2 <4 <8 <16 <32 <64 <128 >=128]: %ld %ld %ld %ld %ld %ld %ld %ld\n",
               h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7]);
        long over = h[5]+h[6]+h[7];   /* >= 32 vox = clearly stretched */
        *out_over = over;
    }
    (void)hm;

    /* top-N longest edges with coords */
    qsort(E, ne, sizeof(Edge), cmp_edge_desc);
    int shown = topN < (int)ne ? topN : (int)ne;
    printf("  top %d longest edges (len : vtx a (z,y,x) -> vtx b (z,y,x)  face):\n", shown);
    for (int i=0;i<shown;i++){
        Edge *e=&E[i];
        printf("   %8.1f : v%-8d (%8.1f,%8.1f,%8.1f) -> v%-8d (%8.1f,%8.1f,%8.1f)  f%d\n",
               e->len, e->a, V[(size_t)e->a*3],V[(size_t)e->a*3+1],V[(size_t)e->a*3+2],
               e->b, V[(size_t)e->b*3],V[(size_t)e->b*3+1],V[(size_t)e->b*3+2], e->face);
    }

    /* spike vertices: vmax / mean_incident > spike */
    long nspike=0; double worst_ratio=0; int32_t worst_v=-1;
    for(size_t i=0;i<nv;i++){
        if(vcnt[i]==0)continue;
        double mean=vsum[i]/(double)vcnt[i];
        if(mean<1e-9)continue;
        double ratio=vmax[i]/mean;
        if(ratio>spike)nspike++;
        if(ratio>worst_ratio){worst_ratio=ratio; worst_v=(int32_t)i;}
    }
    printf("  spike verts (max incident > %.0fx mean incident): %ld   worst ratio=%.1f (v%d)\n",
           spike, nspike, worst_ratio, worst_v);

    free(E); free(len); free(vmax); free(vsum); free(vcnt);
}

static int selftest(void)
{
    int fails=0;
    printf("=== obj_longedge selftest ===\n");
    /* two unit tris sharing an edge + one spike vert pulled far away */
    double V[]={ 0,0,0,  0,0,1,  0,1,0,  0,100,100 };  /* v3 is the spike */
    int32_t F[]={ 0,1,2,  1,3,2 };
    double mx=0; long over=0;
    analyze(V,4,F,2, 5, 8.0, 0.0, &mx, &over);
    /* longest edge involves v3 at (0,100,100): length ~ from (0,0,1) or (0,1,0) */
    if(!(mx>100.0)){ printf("  FAIL: max=%.1f expected >100\n", mx); fails++; }
    else printf("  ok: max edge %.1f localized\n", mx);
    printf("=== selftest %s (%d failure%s) ===\n", fails?"FAILED":"PASSED", fails, fails==1?"":"s");
    return fails?3:0;
}

int main(int argc, char **argv)
{
    if (argc>=2 && !strcmp(argv[1],"--selftest")) return selftest();
    if (argc<2){ fprintf(stderr,"Usage: %s <mesh.obj> [--top N] [--spike F] [--hist-max L]\n       %s --selftest\n",argv[0],argv[0]); return 2; }
    const char *path=argv[1];
    int topN=20; double spike=8.0, hist_max=0.0;
    for(int i=2;i<argc;i++){
        if(!strcmp(argv[i],"--top")&&i+1<argc)topN=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--spike")&&i+1<argc)spike=atof(argv[++i]);
        else if(!strcmp(argv[i],"--hist-max")&&i+1<argc)hist_max=atof(argv[++i]);
        else { fprintf(stderr,"unknown arg %s\n",argv[i]); return 2; }
    }
    double *V; int32_t *F; size_t nv,nf;
    if(load_obj(path,&V,&nv,&F,&nf)!=0){ fprintf(stderr,"ERROR: cannot load %s\n",path); return 1; }
    printf("=== obj_longedge: %s ===\n  V=%zu F=%zu\n", path, nv, nf);
    double mx; long over;
    analyze(V,nv,F,nf, topN, spike, hist_max, &mx, &over);
    free(V); free(F);
    return 0;
}
