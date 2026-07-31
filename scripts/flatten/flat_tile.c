/* flat_tile.c -- split a wide flat unroll ("v u v 0 [r g b]") into N horizontal
 * u-tiles, preserving per-vertex colors, so each tile renders at a readable
 * aspect ratio. A face is emitted into the tile its centroid-u falls in
 * (each face once). Vertices are compacted per tile.
 *
 * Usage:
 *   flat_tile <in.obj> <out_prefix> <ntiles>
 *   flat_tile --selftest
 * Writes <out_prefix>_00.obj .. <out_prefix>_(N-1).obj and prints each tile's
 * u-range + face count.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct { double u, v; float r, g, b; int has_col; } Vtx;
typedef struct { Vtx *p; size_t n, cap; } VArr;
typedef struct { int32_t *p; size_t n, cap; } IArr;
static void vpush(VArr *a, Vtx x){ if(a->n==a->cap){a->cap=a->cap?a->cap*2:1<<16;a->p=realloc(a->p,a->cap*sizeof(Vtx)); if(!a->p){fprintf(stderr,"OOM\n");exit(1);}} a->p[a->n++]=x; }
static void ipush(IArr *a, int32_t x){ if(a->n==a->cap){a->cap=a->cap?a->cap*2:1<<16;a->p=realloc(a->p,a->cap*sizeof(int32_t)); if(!a->p){fprintf(stderr,"OOM\n");exit(1);}} a->p[a->n++]=x; }

static int load(const char *path, VArr *V, IArr *F)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    char line[8192];
    while (fgets(line, sizeof line, f)) {
        if (line[0]=='v' && line[1]==' ') {
            Vtx x; memset(&x,0,sizeof x);
            double a,b,c,r,g,bl;
            int n = sscanf(line+2, "%lf %lf %lf %lf %lf %lf", &a,&b,&c,&r,&g,&bl);
            if (n>=3){ x.u=a; x.v=b; if(n>=6){x.r=(float)r;x.g=(float)g;x.b=(float)bl;x.has_col=1;} vpush(V,x); }
        } else if (line[0]=='f' && line[1]==' ') {
            int idx[64], ni=0; const char *p=line+1;
            while(*p && ni<64){ while(*p==' '||*p=='\t')p++; if(*p=='\0'||*p=='\n'||*p=='\r')break;
                long vi=strtol(p,(char**)&p,10); while(*p&&*p!=' '&&*p!='\t'&&*p!='\n'&&*p!='\r')p++;
                idx[ni++]= vi>0? (int32_t)(vi-1) : (int32_t)((long)V->n+vi); }
            for(int i=1;i+1<ni;i++){ ipush(F,idx[0]); ipush(F,idx[i]); ipush(F,idx[i+1]); }
        }
    }
    fclose(f);
    return 0;
}

static int write_tile(const char *path, const VArr *V, const IArr *F,
                      double u0, double u1, size_t *out_nf)
{
    size_t nv = V->n, nf = F->n/3;
    int32_t *remap = malloc(nv*sizeof(int32_t));
    for (size_t i=0;i<nv;i++) remap[i]=-1;
    FILE *f = fopen(path, "wb");
    if(!f){ free(remap); return -1; }
    /* first pass: which verts are used by in-tile faces */
    size_t nfout=0, nvout=0;
    for (size_t t=0;t<nf;t++){
        int32_t a=F->p[t*3],b=F->p[t*3+1],c=F->p[t*3+2];
        double cu=(V->p[a].u+V->p[b].u+V->p[c].u)/3.0;
        if(cu<u0||cu>=u1) continue;
        int32_t vs[3]={a,b,c};
        for(int k=0;k<3;k++) if(remap[vs[k]]<0) remap[vs[k]]=(int32_t)(nvout++);
        nfout++;
    }
    /* emit verts in original order among used */
    int32_t *order = malloc(nvout*sizeof(int32_t));
    for(size_t i=0;i<nv;i++) if(remap[i]>=0) order[remap[i]]=(int32_t)i;
    for(size_t i=0;i<nvout;i++){ const Vtx *x=&V->p[order[i]];
        if(x->has_col) fprintf(f,"v %.4f %.4f 0 %.3f %.3f %.3f\n",x->u,x->v,x->r,x->g,x->b);
        else fprintf(f,"v %.4f %.4f 0\n",x->u,x->v); }
    for (size_t t=0;t<nf;t++){
        int32_t a=F->p[t*3],b=F->p[t*3+1],c=F->p[t*3+2];
        double cu=(V->p[a].u+V->p[b].u+V->p[c].u)/3.0;
        if(cu<u0||cu>=u1) continue;
        fprintf(f,"f %d %d %d\n",remap[a]+1,remap[b]+1,remap[c]+1);
    }
    fclose(f); free(remap); free(order);
    *out_nf=nfout;
    return 0;
}

static int selftest(void)
{
    /* two tris, one at u~1, one at u~11; 2 tiles [0,10),[10,20) -> 1 face each */
    const char *in="flat_tile_st_in.obj", *pre="flat_tile_st";
    FILE *f=fopen(in,"wb");
    fprintf(f,"v 1 0 0 1 0 0\nv 1 1 0 1 0 0\nv 2 0 0 1 0 0\n");
    fprintf(f,"v 11 0 0 0 1 0\nv 11 1 0 0 1 0\nv 12 0 0 0 1 0\n");
    fprintf(f,"f 1 2 3\nf 4 5 6\n"); fclose(f);
    VArr V={0}; IArr F={0}; load(in,&V,&F);
    char p0[64],p1[64]; snprintf(p0,64,"%s_00.obj",pre); snprintf(p1,64,"%s_01.obj",pre);
    size_t n0,n1; write_tile(p0,&V,&F,0,10,&n0); write_tile(p1,&V,&F,10,20,&n1);
    int ok = (n0==1 && n1==1);
    /* check color preserved in tile 0 */
    FILE *t=fopen(p0,"rb"); char ln[256]; int red=0;
    while(fgets(ln,sizeof ln,t)) if(strstr(ln,"1.000 0.000 0.000")) red=1;
    fclose(t); ok = ok && red;
    remove(in); remove(p0); remove(p1);
    printf("=== flat_tile selftest %s (n0=%zu n1=%zu redcol=%d) ===\n", ok?"PASSED":"FAILED",n0,n1,red);
    free(V.p); free(F.p);
    return ok?0:3;
}

int main(int argc, char **argv)
{
    if (argc>=2 && !strcmp(argv[1],"--selftest")) return selftest();
    if (argc<4){ fprintf(stderr,"Usage: %s <in.obj> <out_prefix> <ntiles>\n       %s --selftest\n",argv[0],argv[0]); return 2; }
    const char *in=argv[1], *pre=argv[2]; int nt=atoi(argv[3]);
    if(nt<1)nt=1;
    VArr V={0}; IArr F={0};
    if(load(in,&V,&F)!=0){ fprintf(stderr,"ERROR: cannot load %s\n",in); return 1; }
    double umin=1e300,umax=-1e300;
    for(size_t i=0;i<V.n;i++){ if(V.p[i].u<umin)umin=V.p[i].u; if(V.p[i].u>umax)umax=V.p[i].u; }
    double span=(umax-umin)/(double)nt;
    printf("=== flat_tile: %s  V=%zu F=%zu  u=[%.0f,%.0f] -> %d tiles ===\n", in, V.n, F.n/3, umin, umax, nt);
    for(int t=0;t<nt;t++){
        double u0=umin+span*t, u1=(t==nt-1)?umax+1.0:umin+span*(t+1);
        char path[1024]; snprintf(path,sizeof path,"%s_%02d.obj",pre,t);
        size_t nf=0; write_tile(path,&V,&F,u0,u1,&nf);
        printf("  tile %02d  u=[%.0f,%.0f]  %zu faces -> %s\n", t, u0, u1, nf, path);
    }
    free(V.p); free(F.p);
    return 0;
}
