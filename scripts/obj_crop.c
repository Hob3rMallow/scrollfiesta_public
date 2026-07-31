/* obj_crop — keep only the faces whose 3 vertices all fall inside a world-space
 * bbox, remap/compact vertices, write a small OBJ. For zooming in on a single
 * BPA tear / seam gap so it is actually visible in a render.
 *
 * Usage:  obj_crop <in.obj> <out.obj> z0 z1 y0 y1 x0 x1
 *         obj_crop --selftest
 *
 * OBJ vertex order is (z,y,x) per repo convention; the bbox args match that.
 * Standalone C (plain malloc) — no repo deps.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { float z, y, x; } V3;

static int inb(V3 v, const double *b)
{
    return v.z>=b[0]&&v.z<=b[1]&&v.y>=b[2]&&v.y<=b[3]&&v.x>=b[4]&&v.x<=b[5];
}

/* Crop: returns #faces kept; writes out.obj. Arrays 1-based faces (OBJ). */
static int crop(const char *in, const char *out, const double *bb,
                long *out_nv, long *out_nf)
{
    FILE *f = fopen(in, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", in); return -1; }
    long vcap = 1<<16, nv = 0;
    V3 *V = (V3 *)malloc((size_t)vcap*sizeof(V3));
    long fcap = 1<<16, nf = 0;
    int *F = (int *)malloc((size_t)fcap*3*sizeof(int));
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (line[0]=='v'&&line[1]==' ') {
            V3 v;
            if (sscanf(line+2, "%f %f %f", &v.z,&v.y,&v.x)!=3) continue;
            if (nv>=vcap){vcap*=2;V=(V3*)realloc(V,(size_t)vcap*sizeof(V3));}
            V[nv++]=v;
        } else if (line[0]=='f'&&line[1]==' ') {
            int a,b,c;
            /* handles "f a b c" and "f a/.. b/.. c/.." */
            if (sscanf(line+2, "%d",&a)!=1) continue;
            char *p=line+2; int idx[3],k=0;
            while (*p&&k<3){ while(*p==' ')p++; if(!*p)break;
                idx[k++]=atoi(p); while(*p&&*p!=' ')p++; }
            if (k<3) continue;
            a=idx[0];b=idx[1];c=idx[2];
            if (nf>=fcap){fcap*=2;F=(int*)realloc(F,(size_t)fcap*3*sizeof(int));}
            F[nf*3+0]=a;F[nf*3+1]=b;F[nf*3+2]=c;nf++;
        }
    }
    fclose(f);

    int *remap = (int *)malloc((size_t)(nv+1)*sizeof(int));
    for (long i=0;i<=nv;i++) remap[i]=0;
    FILE *o = fopen(out, "w");
    if (!o) { fprintf(stderr,"cannot write %s\n",out); return -1; }
    /* pass 1: mark used verts of kept faces */
    long kf=0;
    for (long i=0;i<nf;i++){
        int a=F[i*3+0],b=F[i*3+1],c=F[i*3+2];
        if(a<1||a>nv||b<1||b>nv||c<1||c>nv)continue;
        if(inb(V[a-1],bb)&&inb(V[b-1],bb)&&inb(V[c-1],bb)){
            remap[a]=remap[b]=remap[c]=1; kf++;
        }
    }
    /* assign compact indices, emit verts */
    long knv=0;
    for (long i=1;i<=nv;i++) if(remap[i]){ knv++; remap[i]=(int)knv;
        fprintf(o,"v %.4f %.4f %.4f\n",(double)V[i-1].z,(double)V[i-1].y,(double)V[i-1].x); }
    /* emit faces */
    for (long i=0;i<nf;i++){
        int a=F[i*3+0],b=F[i*3+1],c=F[i*3+2];
        if(a<1||a>nv||b<1||b>nv||c<1||c>nv)continue;
        if(inb(V[a-1],bb)&&inb(V[b-1],bb)&&inb(V[c-1],bb))
            fprintf(o,"f %d %d %d\n",remap[a],remap[b],remap[c]);
    }
    fclose(o);
    free(V);free(F);free(remap);
    *out_nv=knv; *out_nf=kf;
    return 0;
}

static int selftest(void)
{
    /* Two triangles; one inside unit box, one far outside. Crop keeps 1. */
    const char *in="_crop_st_in.obj",*out="_crop_st_out.obj";
    FILE *f=fopen(in,"w");
    fprintf(f,"v 0 0 0\nv 0 0 1\nv 0 1 0\n");      /* tri A inside */
    fprintf(f,"v 100 100 100\nv 100 100 101\nv 100 101 100\n"); /* tri B outside */
    fprintf(f,"f 1 2 3\nf 4 5 6\n");
    fclose(f);
    double bb[6]={-1,2,-1,2,-1,2};
    long nv=0,nf=0;
    if(crop(in,out,bb,&nv,&nf)!=0){printf("SELFTEST FAIL: crop error\n");return 1;}
    remove(in);remove(out);
    if(nv!=3||nf!=1){printf("SELFTEST FAIL: kept nv=%ld nf=%ld (want 3,1)\n",nv,nf);return 1;}
    printf("SELFTEST OK (kept 1 face / 3 verts)\n");
    return 0;
}

int main(int argc, char **argv)
{
    if (argc>=2 && !strcmp(argv[1],"--selftest")) return selftest();
    if (argc<9){ fprintf(stderr,"usage: %s <in.obj> <out.obj> z0 z1 y0 y1 x0 x1\n"
                                 "   or: %s --selftest\n",argv[0],argv[0]); return 1; }
    double bb[6]; for(int k=0;k<6;k++) bb[k]=atof(argv[3+k]);
    long nv=0,nf=0;
    if (crop(argv[1],argv[2],bb,&nv,&nf)!=0) return 1;
    printf("cropped %s -> %s : %ld verts, %ld faces  (bbox z[%.0f,%.0f] y[%.0f,%.0f] x[%.0f,%.0f])\n",
           argv[1],argv[2],nv,nf,bb[0],bb[1],bb[2],bb[3],bb[4],bb[5]);
    return 0;
}
