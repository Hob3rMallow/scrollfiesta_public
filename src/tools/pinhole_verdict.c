/* pinhole_verdict.c — explain, per single-triangle boundary pinhole in an OBJ,
 * why ball_pivot.c::fill_pinholes would or would not have closed it.
 *
 * For each boundary loop of length <= maxloop (default 3 = a one-triangle hole)
 * it replays the two gates that actually decide the fill:
 *
 *   1. TRACER gate (fill_pinholes): the loop is only a candidate if it traces as
 *      a simple cycle, i.e. every boundary vertex has boundary-degree exactly 2.
 *      A pinch / shared-boundary vertex (degree != 2) makes simple=0 and the hole
 *      is skipped regardless of geometry.
 *
 *   2. GEOMETRY gate (pinhole_fillable): on the surface-consistent ball side the
 *      rho-ball must be empty; a borderline-cocircular intruder is judged by the
 *      exact translated-double in-circle test (robust_inside_circum, copied from
 *      ball_pivot.c verbatim).
 *
 * The dump carries no MLS normals, so the surface normal is taken from the face
 * WINDING (BPA winds every triangle so its normal follows +N), which is the same
 * side the normal gate uses. For each hole the tool prints both ball sides, the
 * nearest intruder + depth, the intruder's connectivity to the hole, and a final
 * verdict naming the responsible gate.
 *
 * Usage:  pinhole_verdict <mesh.obj> [rho=1.5] [maxloop=3]
 *
 * Standalone C99; no project deps so it builds as its own tiny vcxproj.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* ---- tiny vec3 helpers (double) ---- */
static void  v_sub(double *o,const double *a,const double *b){o[0]=a[0]-b[0];o[1]=a[1]-b[1];o[2]=a[2]-b[2];}
static void  v_cross(double *o,const double *a,const double *b){o[0]=a[1]*b[2]-a[2]*b[1];o[1]=a[2]*b[0]-a[0]*b[2];o[2]=a[0]*b[1]-a[1]*b[0];}
static double v_dot(const double *a,const double *b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static double v_norm(const double *a){return sqrt(v_dot(a,a));}

/* ---- exact coplanar in-circle test (mirrors the FIXED ball_pivot.c) ----
 * Signed value: >0 == query P strictly inside the circumcircle of A,B,C.
 * Uses an ORTHONORMAL in-plane basis (u,v) from the plane normal so the lifted-
 * paraboloid determinant tests the true circle, not the sheared ellipse an axis-
 * aligned coordinate drop would produce on a tilted plane. Winding-normalised. */
static double robust_incircle_value(const double *a,const double *b,const double *c,
                                    const double *p,const double *nrm)
{
    double n[3]={nrm[0],nrm[1],nrm[2]};
    double nl=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    if(nl<1e-20){ double e1[3],e2[3]; v_sub(e1,b,a); v_sub(e2,c,a); v_cross(n,e1,e2);
        nl=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); if(nl<1e-20) return 0.0; }
    n[0]/=nl;n[1]/=nl;n[2]/=nl;
    double t[3]={1,0,0}; if(fabs(n[0])>0.9){ t[0]=0; t[1]=1; }   /* keep t away from n */
    double u[3],v[3]; v_cross(u,t,n); double ul=v_norm(u); u[0]/=ul;u[1]/=ul;u[2]/=ul; v_cross(v,n,u);
    double db[3],dc[3],dp[3]; v_sub(db,b,a); v_sub(dc,c,a); v_sub(dp,p,a);
    double bx=v_dot(db,u),by=v_dot(db,v), cx=v_dot(dc,u),cy=v_dot(dc,v), px=v_dot(dp,u),py=v_dot(dp,v);
    double orient=bx*cy-by*cx;
    double bl=bx*bx+by*by, cl=cx*cx+cy*cy, pl=px*px+py*py;
    double det=bx*(cy*pl-cl*py)-by*(cx*pl-cl*px)+bl*(cx*py-cy*px);
    double val=(orient>=0.0)?det:-det;     /* <0 == inside (lifted paraboloid) */
    return -val;                            /* flip so >0 == inside */
}
static int robust_inside_circum(const double *a,const double *b,const double *c,
                                const double *p,const double *nrm)
{ return robust_incircle_value(a,b,c,p,nrm) > 0.0; }

/* ---- growable arrays ---- */
typedef struct { double *p; size_t n, cap; } DArr;   /* xyz triples flattened */
typedef struct { int    *p; size_t n, cap; } IArr;
static void darr_push3(DArr *d,double x,double y,double z){
    if(d->n+3>d->cap){ d->cap=d->cap?d->cap*2:3072; d->p=realloc(d->p,d->cap*sizeof(double)); }
    d->p[d->n++]=x; d->p[d->n++]=y; d->p[d->n++]=z;
}
static void iarr_push3(IArr *a,int x,int y,int z){
    if(a->n+3>a->cap){ a->cap=a->cap?a->cap*2:3072; a->p=realloc(a->p,a->cap*sizeof(int)); }
    a->p[a->n++]=x; a->p[a->n++]=y; a->p[a->n++]=z;
}

/* ---- undirected-edge hash map: key=(lo<<32|hi) -> incident face count ---- */
typedef struct { uint64_t *key; int *val; size_t cap, n; } EdgeMap;
static void em_init(EdgeMap *m,size_t hint){ m->cap=1; while(m->cap<hint*2) m->cap<<=1;
    m->key=calloc(m->cap,sizeof(uint64_t)); m->val=calloc(m->cap,sizeof(int)); m->n=0; }
static int *em_slot(EdgeMap *m,int a,int b){
    uint64_t lo=(a<b)?(uint64_t)a:(uint64_t)b, hi=(a<b)?(uint64_t)b:(uint64_t)a;
    uint64_t k=((lo+1)<<32)|(hi+1);                 /* +1 so key 0 = empty */
    size_t h=(size_t)((k*1099511628211ULL)&(m->cap-1));
    for(;;){ if(m->key[h]==0){ m->key[h]=k; m->n++; return &m->val[h]; }
             if(m->key[h]==k) return &m->val[h];
             h=(h+1)&(m->cap-1); }
}
static int em_get(EdgeMap *m,int a,int b){
    uint64_t lo=(a<b)?(uint64_t)a:(uint64_t)b, hi=(a<b)?(uint64_t)b:(uint64_t)a;
    uint64_t k=((lo+1)<<32)|(hi+1);
    size_t h=(size_t)((k*1099511628211ULL)&(m->cap-1));
    for(;;){ if(m->key[h]==0) return 0; if(m->key[h]==k) return m->val[h]; h=(h+1)&(m->cap-1); }
}

/* ---- self-test: is robust_inside_circum's sign correct? --------------------
 * Ground truth comes from the unambiguous geometric test: project the query
 * point onto ABC's plane and compare its distance from the circumcenter to the
 * circumradius. We feed several triangles (both windings, axis-aligned and
 * tilted planes) and assert robust_inside_circum agrees. */
static int circum3d(const double*a,const double*b,const double*c,double*cc_out,double*rc_out,double*nu_out){
    double ap[3],bp[3]; v_sub(ap,a,c); v_sub(bp,b,c);
    double m[3]; v_cross(m,ap,bp); double m2=v_dot(m,m); if(m2<1e-20) return 0;
    double la=v_dot(ap,ap), lb=v_dot(bp,bp);
    double lhs[3]={la*bp[0]-lb*ap[0], la*bp[1]-lb*ap[1], la*bp[2]-lb*ap[2]};
    double pp[3]; v_cross(pp,lhs,m); pp[0]/=2*m2; pp[1]/=2*m2; pp[2]/=2*m2;
    cc_out[0]=c[0]+pp[0]; cc_out[1]=c[1]+pp[1]; cc_out[2]=c[2]+pp[2];
    *rc_out=v_norm(pp); double ml=sqrt(m2); nu_out[0]=m[0]/ml; nu_out[1]=m[1]/ml; nu_out[2]=m[2]/ml;
    return 1;
}
static int gt_inside(const double*a,const double*b,const double*c,const double*p){
    double cc[3],rc,nu[3]; if(!circum3d(a,b,c,cc,&rc,nu)) return -1;
    double pc[3]; v_sub(pc,p,cc); double hgt=v_dot(pc,nu);
    double inplane=sqrt(fmax(0.0, v_dot(pc,pc)-hgt*hgt));
    if(fabs(inplane-rc)<1e-6) return -1;        /* on the circle: skip (tie) */
    return inplane<rc ? 1 : 0;
}
static int selftest(void){
    /* case list: {A,B,C, P, normal} */
    double cases[][15] = {
        /* CCW xy triangle, P=circumcenter (inside) */
        {0,0,0, 1,0,0, 0,1,0,  0.5,0.5,0,   0,0,1},
        /* same, P far outside */
        {0,0,0, 1,0,0, 0,1,0,  3,3,0,       0,0,1},
        /* CW winding (swap B,C), P=center (inside) */
        {0,0,0, 0,1,0, 1,0,0,  0.5,0.5,0,   0,0,1},
        /* tilted plane (normal ~ (1,1,1)), P near centroid (inside) */
        {0,0,0, 1,0,0, 0,1,1,  0.3,0.3,0.3, 0,0,0},
        /* tilted plane, P pushed far along an in-plane direction (outside) */
        {0,0,0, 1,0,0, 0,1,1,  5,-3,2,      0,0,0},
    };
    int n=(int)(sizeof(cases)/sizeof(cases[0]));
    int fails=0;
    for(int i=0;i<n;i++){
        double *a=&cases[i][0], *b=&cases[i][3], *c=&cases[i][6], *p=&cases[i][9], *nrm=&cases[i][12];
        double nbuf[3]={nrm[0],nrm[1],nrm[2]};
        if(v_dot(nbuf,nbuf)<1e-12){ double e1[3],e2[3]; v_sub(e1,b,a); v_sub(e2,c,a); v_cross(nbuf,e1,e2); }
        int gt=gt_inside(a,b,c,p);
        int rv=robust_inside_circum(a,b,c,p,nbuf);
        double det=robust_incircle_value(a,b,c,p,nbuf);
        const char *ok = (gt<0)?"tie/skip":(gt==rv?"ok":"MISMATCH");
        printf("[selftest %d] ground-truth=%s  robust=%s  det=%+.3e  -> %s\n",
               i, gt<0?"tie":(gt?"INSIDE":"OUTSIDE"), rv?"INSIDE":"OUTSIDE", det, ok);
        if(gt>=0 && gt!=rv) fails++;
    }
    printf("[selftest] %s (%d mismatch%s)\n", fails?"FAILED — robust_inside_circum sign is wrong":"PASSED",
           fails, fails==1?"":"es");
    return fails?1:0;
}

/* ---- pin band (mirrors src/common/pin_mask.c::Pin_recompute_one) ----
 * A vertex is pinned iff, in CUBE-LOCAL coords (world - cube_origin), any axis
 * is within PIN_EPS of a cube face (0 or CUBE_D). These verts hold their raw MC
 * positions for cross-cube welding; the holefill loop-skip leaves any loop that
 * touches one (the assumption: a pinned vert is a seam edge the weld will close).
 * Constants must match pipeline_constants.h (HALO_PIN_EPS) and main.c (cube_D). */
#define PIN_EPS  1.5
#define CUBE_D   128.0
static const char *FACE_NAME[6] = {"z-","z+","y-","y+","x-","x+"};

/* Parse cube origin (z,y,x) from a "...z%d_y%d_x%d..." path. Mirrors
 * dump_obj.c::parse_cube_origin_zyx but scans anywhere in the basename. */
static int parse_origin_from_path(const char *path,double org[3]){
    org[0]=org[1]=org[2]=0.0;
    const char *b=path, *s;
    for(s=path; *s; s++) if(*s=='/'||*s=='\\') b=s+1;   /* basename */
    const char *z=strstr(b,"z");
    for(; z && *z; z=strstr(z+1,"z")){
        int iz=0,iy=0,ix=0;
        if(sscanf(z,"z%d_y%d_x%d",&iz,&iy,&ix)==3){
            org[0]=iz; org[1]=iy; org[2]=ix; return 0;
        }
    }
    return -1;
}
/* face bitmask (6 bits, see FACE_NAME) for a world vert given the cube origin */
static int pin_faces(const double *wv,const double *org){
    double l[3]={wv[0]-org[0],wv[1]-org[1],wv[2]-org[2]};
    int m=0;
    for(int ax=0;ax<3;ax++){
        if(l[ax] < PIN_EPS)          m |= 1<<(ax*2);     /* low  face  */
        if(l[ax] > CUBE_D - PIN_EPS) m |= 1<<(ax*2+1);   /* high face  */
    }
    return m;
}

int main(int argc,char**argv)
{
    if(argc>=2 && strcmp(argv[1],"--selftest")==0) return selftest();
    if(argc<2){ fprintf(stderr,"Usage: %s <mesh.obj> [rho=1.5] [maxloop=3]\n"
                                "       %s <mesh.obj> --pins [--seam x+,z-]   (write bright-red pin OBJ)\n"
                                "       %s --selftest\n",argv[0],argv[0],argv[0]); return 2; }
    const char *path=argv[1];
    double rho = (argc>2 && argv[2][0]!='-')? atof(argv[2]) : 1.5;
    int maxloop = (argc>3 && argv[3][0]!='-')? atoi(argv[3]) : 3;
    double degen_thresh = -1.0;   /* --degen <deg>: list triangles below <deg> min-angle */
    int pinch_report = 0;         /* --pinch: characterize boundary pinch points (bdeg>2) */
    int pins_report = 0;          /* --pins: pin-band analysis + bright-red colored OBJ */
    const char *seam_arg = NULL;  /* --seam <faces>: which faces are real seams (e.g. x+,z-) */
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i],"--degen")==0 && i+1<argc) degen_thresh=atof(argv[i+1]);
        if(strcmp(argv[i],"--pinch")==0) pinch_report=1;
        if(strcmp(argv[i],"--pins")==0)  pins_report=1;
        if(strcmp(argv[i],"--seam")==0 && i+1<argc) seam_arg=argv[i+1];
    }

    FILE *f=fopen(path,"r");
    if(!f){ fprintf(stderr,"ERROR: cannot open %s\n",path); return 1; }

    DArr V={0}; IArr Fc={0};
    char line[4096];
    while(fgets(line,sizeof(line),f)){
        if(line[0]=='v'&&line[1]==' '){
            double x=0,y=0,z=0;
            if(sscanf(line+2,"%lf %lf %lf",&x,&y,&z)==3) darr_push3(&V,x,y,z);
        } else if(line[0]=='f'&&line[1]==' '){
            int idx[16]; int ni=0;
            const char *p=line+1;
            while(*p && ni<16){
                while(*p==' '||*p=='\t') p++;
                if(*p=='\0'||*p=='\n'||*p=='\r') break;
                int vi=(int)strtol(p,(char**)&p,10);
                while(*p && *p!=' ' && *p!='\t' && *p!='\n' && *p!='\r') p++; /* skip /vt/vn */
                if(vi<0) vi=(int)(V.n/3)+vi; else vi-=1;
                idx[ni++]=vi;
            }
            for(int i=1;i+1<ni;i++) iarr_push3(&Fc,idx[0],idx[i],idx[i+1]);
        }
    }
    fclose(f);
    size_t nv=V.n/3, nf=Fc.n/3;
    printf("=== pinhole_verdict: %s ===\n", path);
    printf("    V=%zu  F=%zu  rho=%.3f  maxloop=%d\n", nv, nf, rho, maxloop);
    if(nv<3||nf<1){ printf("    (empty mesh)\n"); return 0; }

    /* per-vertex incident faces (CSR) */
    int *vf_cnt=calloc(nv,sizeof(int));
    for(size_t i=0;i<Fc.n;i++) vf_cnt[Fc.p[i]]++;
    size_t *vf_off=malloc((nv+1)*sizeof(size_t)); vf_off[0]=0;
    for(size_t i=0;i<nv;i++) vf_off[i+1]=vf_off[i]+(size_t)vf_cnt[i];
    int *vf_idx=malloc(vf_off[nv]*sizeof(int));
    { size_t *cur=malloc(nv*sizeof(size_t)); memcpy(cur,vf_off,nv*sizeof(size_t));
      for(size_t fi=0;fi<nf;fi++) for(int k=0;k<3;k++){ int v=Fc.p[fi*3+k]; vf_idx[cur[v]++]=(int)fi; }
      free(cur); }

    /* undirected edge incidence counts */
    EdgeMap em; em_init(&em, nf*3);
    for(size_t fi=0;fi<nf;fi++){ int a=Fc.p[fi*3],b=Fc.p[fi*3+1],c=Fc.p[fi*3+2];
        (*em_slot(&em,a,b))++; (*em_slot(&em,b,c))++; (*em_slot(&em,c,a))++; }

    /* boundary degree per vertex + directed next[] for loop tracing */
    int *bdeg=calloc(nv,sizeof(int));
    int *next=malloc(nv*sizeof(int)); for(size_t i=0;i<nv;i++) next[i]=-1;
    for(size_t fi=0;fi<nf;fi++){ int v[3]={Fc.p[fi*3],Fc.p[fi*3+1],Fc.p[fi*3+2]};
        for(int e=0;e<3;e++){ int u=v[e], w=v[(e+1)%3];
            if(em_get(&em,u,w)==1){ next[u]=w; } } }
    /* boundary degree = number of undirected boundary edges touching the vertex.
     * Iterate the edge map so BOTH endpoints of each boundary edge are counted
     * (a simple-loop vertex therefore has bdeg==2, matching fill_pinholes' deg[]
     * which does deg[ba]++; deg[bb]++ per boundary edge). */
    for(size_t i=0;i<em.cap;i++){
        if(em.key[i] && em.val[i]==1){
            int lo=(int)((em.key[i]>>32)-1), hi=(int)((em.key[i]&0xffffffffULL)-1);
            bdeg[lo]++; bdeg[hi]++;
        }
    }

    /* --pinch: boundary verts with bdeg>2 are pinch points (two boundary loops
     * meeting at one vertex — what HoleFill's CDT rejects as "duplicate vert").
     * Histogram their distance to the nearest cube-face plane (faces sit at
     * multiples of 128 vox in world coords) to tell chewed-at-seam from general
     * BPA raggedness. */
    if(pinch_report){
        long np=0; long fh[6]={0,0,0,0,0,0}; /* <1 1-3 3-6 6-12 12-32 >32 vox */
        double nearsum=0;
        for(size_t i=0;i<nv;i++){
            if(bdeg[i]<=2) continue;
            np++;
            double fd=1e30;
            for(int k=0;k<3;k++){ double c=V.p[i*3+k]; double m=c-floor(c/128.0)*128.0;
                double d=m<128.0-m?m:128.0-m; if(d<fd) fd=d; }
            nearsum+=fd;
            int b=(fd<1)?0:(fd<3)?1:(fd<6)?2:(fd<12)?3:(fd<32)?4:5; fh[b]++;
        }
        printf("  pinch: %ld boundary verts with bdeg>2 (pinch points). dist-to-nearest-128-plane:\n"
               "         [<1]=%ld [1-3]=%ld [3-6]=%ld [6-12]=%ld [12-32]=%ld [>32]=%ld  mean=%.2f vox\n",
               np, fh[0],fh[1],fh[2],fh[3],fh[4],fh[5], np?nearsum/(double)np:0.0);
    }

    /* --pins: replay the holefill pin-skip gate. Recompute the pin band from
     * cube geometry (pin_mask.c), report which faces are pinned and — per
     * boundary loop — how many verts are pinned and on which faces, then write
     * a bright-red pin-colored OBJ. A loop pinned ONLY on a face that is NOT a
     * real seam (--seam) is a hole the weld will never close: wrongly skipped. */
    if(pins_report){
        double org[3];
        if(parse_origin_from_path(path,org)!=0){
            printf("  pins: ERROR — could not parse cube origin (z%%d_y%%d_x%%d) from path; "
                   "cannot recompute the pin band.\n");
        } else {
            int seam_mask=0;
            if(seam_arg) for(int fa=0;fa<6;fa++) if(strstr(seam_arg,FACE_NAME[fa])) seam_mask|=1<<fa;

            /* per-vertex pin face mask */
            int *pmask=calloc(nv,sizeof(int));
            size_t npin=0; long face_cnt[6]={0,0,0,0,0,0};
            for(size_t i=0;i<nv;i++){
                int m=pin_faces(&V.p[i*3],org); pmask[i]=m;
                if(m){ npin++; for(int fa=0;fa<6;fa++) if(m&(1<<fa)) face_cnt[fa]++; }
            }
            printf("  pins: cube origin z=%.0f y=%.0f x=%.0f  eps=%.1f D=%.0f%s%s\n",
                   org[0],org[1],org[2],PIN_EPS,CUBE_D,
                   seam_arg?"  seams=":"", seam_arg?seam_arg:"");
            printf("        pinned %zu/%zu verts (%.1f%%)  per-face[z- z+ y- y+ x- x+]= %ld %ld %ld %ld %ld %ld\n",
                   npin,nv,100.0*(double)npin/(double)nv,
                   face_cnt[0],face_cnt[1],face_cnt[2],face_cnt[3],face_cnt[4],face_cnt[5]);

            /* per boundary loop: walk next[] (no length cap), tally pins/faces */
            uint8_t *pv=calloc(nv,1);
            printf("        boundary loops (pinned verts / length, by face):\n");
            int li=0;
            for(size_t s=0;s<nv;s++){
                if(next[s]<0||pv[s]) continue;
                int cur=(int)s; long len=0,pin=0; int lmask=0,pinmask=0;
                while(cur>=0 && !pv[cur]){ pv[cur]=1; len++;
                    if(pmask[cur]){ pin++; pinmask|=pmask[cur]; } lmask|=pmask[cur];
                    cur=next[cur]; }
                (void)lmask;
                char fb[32]={0}; for(int fa=0;fa<6;fa++) if(pinmask&(1<<fa)){ strcat(fb,FACE_NAME[fa]); strcat(fb," "); }
                int free_only = (seam_arg && pin>0 && (pinmask & ~seam_mask) && !(pinmask & seam_mask));
                printf("          loop %-2d  len=%-4ld pinned=%-4ld %s%s%s%s\n",
                       li, len, pin, pin?"faces=[ ":"(no pins)", pin?fb:"", pin?"]":"",
                       free_only?"  <-- pinned only on FREE (non-seam) faces: weld won't close it":"");
                li++;
            }
            free(pv);

            /* write bright-red pin-colored OBJ next to the input */
            char outp[1200]; snprintf(outp,sizeof(outp),"%.1190s.pins.obj",path);
            FILE *of=fopen(outp,"w");
            if(of){
                for(size_t i=0;i<nv;i++){
                    double r = pmask[i]?1.0:0.55, g = pmask[i]?0.0:0.55, bcol = pmask[i]?0.0:0.55;
                    fprintf(of,"v %.6f %.6f %.6f %.4f %.4f %.4f\n",
                            V.p[i*3],V.p[i*3+1],V.p[i*3+2], r,g,bcol);
                }
                for(size_t fi=0;fi<nf;fi++)
                    fprintf(of,"f %d %d %d\n",Fc.p[fi*3]+1,Fc.p[fi*3+1]+1,Fc.p[fi*3+2]+1);
                fclose(of);
                printf("        wrote pin-colored OBJ -> %s   (pinned=bright red, else gray)\n",outp);
            } else {
                printf("        ERROR: could not write %s\n",outp);
            }
            free(pmask);
        }
    }

    /* ---- topology summary (components via union-find, non-manifold, genus) ---- */
    int *uf=malloc(nv*sizeof(int)); for(size_t i=0;i<nv;i++) uf[i]=(int)i;
    for(size_t fi=0;fi<nf;fi++){
        int vv[3]={Fc.p[fi*3],Fc.p[fi*3+1],Fc.p[fi*3+2]};
        int r0=vv[0]; while(uf[r0]!=r0){uf[r0]=uf[uf[r0]];r0=uf[r0];}
        for(int k=1;k<3;k++){ int r=vv[k]; while(uf[r]!=r){uf[r]=uf[uf[r]];r=uf[r];} if(r!=r0) uf[r]=r0; }
    }
    size_t refV=0, comps=0, nmE=0, bE=0;
    for(size_t i=0;i<em.cap;i++) if(em.key[i]){ if(em.val[i]==1) bE++; else if(em.val[i]>2) nmE++; }
    for(size_t i=0;i<nv;i++) if(vf_cnt[i]>0){ refV++; int r=(int)i; while(uf[r]!=r) r=uf[r]; if(r==(int)i) comps++; }

    /* ---- triangle quality: min interior angle per face (slivers = over-dense BPA) ---- */
    long qhist[7]={0,0,0,0,0,0,0};   /* min-angle buckets: <5 5-10 10-20 20-30 30-40 40-50 50-60 */
    double minang_sum=0.0, elen_sum=0.0; long nbad10=0, nbad30=0;
    for(size_t fi=0;fi<nf;fi++){
        double *p0=&V.p[Fc.p[fi*3]*3], *p1=&V.p[Fc.p[fi*3+1]*3], *p2=&V.p[Fc.p[fi*3+2]*3];
        double e0[3],e1[3],e2[3]; v_sub(e0,p1,p0); v_sub(e1,p2,p1); v_sub(e2,p0,p2);
        double a=v_norm(e0), b=v_norm(e1), c=v_norm(e2);
        elen_sum += a+b+c;
        double mn;
        if(a>1e-12 && b>1e-12 && c>1e-12){
            double ca=(b*b+c*c-a*a)/(2*b*c); if(ca<-1)ca=-1; if(ca>1)ca=1;
            double cb=(a*a+c*c-b*b)/(2*a*c); if(cb<-1)cb=-1; if(cb>1)cb=1;
            double cc=(a*a+b*b-c*c)/(2*a*b); if(cc<-1)cc=-1; if(cc>1)cc=1;
            double A=acos(ca)*57.2957795, B=acos(cb)*57.2957795, C=acos(cc)*57.2957795;
            mn=A; if(B<mn)mn=B; if(C<mn)mn=C;
        } else mn=0.0;
        minang_sum += mn;
        if(mn<10.0) nbad10++; if(mn<30.0) nbad30++;
        int bk=(mn<5)?0:(mn<10)?1:(mn<20)?2:(mn<30)?3:(mn<40)?4:(mn<50)?5:6;
        qhist[bk]++;
        if(degen_thresh>=0.0 && mn<degen_thresh){
            printf("  DEGEN f%zu min-angle=%.2f deg edges=(%.3f,%.3f,%.3f) area2=%.4f\n"
                   "        v%d=(%.2f,%.2f,%.2f) v%d=(%.2f,%.2f,%.2f) v%d=(%.2f,%.2f,%.2f)\n",
                   fi, mn, a,b,c, /* twice area via cross */
                   (double)(v_norm((double[3]){ (p1[1]-p0[1])*(p2[2]-p0[2])-(p1[2]-p0[2])*(p2[1]-p0[1]),
                                                (p1[2]-p0[2])*(p2[0]-p0[0])-(p1[0]-p0[0])*(p2[2]-p0[2]),
                                                (p1[0]-p0[0])*(p2[1]-p0[1])-(p1[1]-p0[1])*(p2[0]-p0[0]) })),
                   Fc.p[fi*3],p0[0],p0[1],p0[2], Fc.p[fi*3+1],p1[0],p1[1],p1[2], Fc.p[fi*3+2],p2[0],p2[1],p2[2]);
        }
    }

    /* trace loops */
    uint8_t *vis=calloc(nv,1);
    int n3=0, would_fill=0, blocked=0, tracer_skip=0, nloops=0;
    for(size_t s=0;s<nv;s++){
        if(next[s]<0||vis[s]) continue;
        int seq[64]; int len=0; int cur=(int)s;
        while(!vis[cur] && next[cur]>=0){ vis[cur]=1; if(len<64) seq[len]=cur; len++; cur=next[cur]; }
        nloops++;
        if(len>maxloop || len<3) continue;
        if(len!=3) continue;                 /* this tool focuses on single-triangle holes */
        n3++;
        int A=seq[0],B=seq[1],C=seq[2];
        double *a=&V.p[A*3], *b=&V.p[B*3], *c=&V.p[C*3];

        /* circumcenter relative to C, circumradius rc */
        double ap[3],bp[3]; v_sub(ap,a,c); v_sub(bp,b,c);
        double m[3]; v_cross(m,ap,bp); double m2=v_dot(m,m), mlen=sqrt(m2);
        printf("\n  --- 3-loop #%d: A=%d B=%d C=%d ---\n",n3,A,B,C);
        printf("      A=(%.2f,%.2f,%.2f) B=(%.2f,%.2f,%.2f) C=(%.2f,%.2f,%.2f)\n",
               a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2]);
        if(m2<1e-20){ printf("      DEGENERATE (collinear) -> no fill\n"); continue; }
        double la=v_dot(ap,ap), lb=v_dot(bp,bp);
        double lhs[3]={la*bp[0]-lb*ap[0], la*bp[1]-lb*ap[1], la*bp[2]-lb*ap[2]};
        double pp[3]; v_cross(pp,lhs,m); pp[0]/=2*m2; pp[1]/=2*m2; pp[2]/=2*m2;
        double cc[3]={c[0]+pp[0],c[1]+pp[1],c[2]+pp[2]}; double rc=v_norm(pp);
        double nu[3]={m[0]/mlen,m[1]/mlen,m[2]/mlen};
        double ctr[3]={(a[0]+b[0]+c[0])/3,(a[1]+b[1]+c[1])/3,(a[2]+b[2]+c[2])/3};

        /* winding-consistent surface normal: mean unit face normal around A,B,C */
        double ns[3]={0,0,0};
        int trio[3]={A,B,C};
        for(int t=0;t<3;t++){ int vtx=trio[t];
            for(size_t j=vf_off[vtx]; j<vf_off[vtx+1]; j++){ int fi=vf_idx[j];
                double *f0=&V.p[Fc.p[fi*3]*3],*f1=&V.p[Fc.p[fi*3+1]*3],*f2=&V.p[Fc.p[fi*3+2]*3];
                double e1[3],e2[3],fn[3]; v_sub(e1,f1,f0); v_sub(e2,f2,f0); v_cross(fn,e1,e2);
                double L=v_norm(fn); if(L>1e-12){ ns[0]+=fn[0]/L; ns[1]+=fn[1]/L; ns[2]+=fn[2]/L; } } }
        double nsl=v_norm(ns); if(nsl>1e-12){ ns[0]/=nsl; ns[1]/=nsl; ns[2]/=nsl; }

        double edgemax=0;
        { double e0[3],e1[3],e2[3]; v_sub(e0,a,b); v_sub(e1,b,c); v_sub(e2,c,a);
          double m0=v_norm(e0),m1=v_norm(e1),m2e=v_norm(e2);
          edgemax=m0>m1?(m0>m2e?m0:m2e):(m1>m2e?m1:m2e); }
        printf("      rc=%.4f rho=%.3f longest-edge=%.4f  bdeg(A,B,C)=(%d,%d,%d)  surf-normal=(%.3f,%.3f,%.3f)\n",
               rc,rho,edgemax,bdeg[A],bdeg[B],bdeg[C],ns[0],ns[1],ns[2]);

        int simple = (bdeg[A]==2 && bdeg[B]==2 && bdeg[C]==2);
        if(rho<rc){ printf("      VERDICT: rho<rc -> no ball touches all three (correctly skipped)\n"); blocked++; continue; }
        double h=sqrt(rho*rho-rc*rc);

        /* both ball sides (context) + identify the surface-consistent side */
        int surf_sgn=0; double surf_O[3]={0,0,0};
        for(int sgn=-1; sgn<=1; sgn+=2){
            double O[3]={cc[0]+sgn*h*nu[0],cc[1]+sgn*h*nu[1],cc[2]+sgn*h*nu[2]};
            double oc[3]; v_sub(oc,O,ctr); int surfside=(v_dot(oc,ns)>0.0);
            int bi=-1; double bd=1e300;
            for(size_t i=0;i<nv;i++){ if((int)i==A||(int)i==B||(int)i==C) continue;
                double dx=V.p[i*3]-O[0],dy=V.p[i*3+1]-O[1],dz=V.p[i*3+2]-O[2]; double dd=dx*dx+dy*dy+dz*dz;
                if(dd<bd){bd=dd;bi=(int)i;} }
            double dist=sqrt(bd);
            printf("      side %c [%s]: nearest other v%d @ %.4f  %s\n", sgn<0?'-':'+',
                   surfside?"SURFACE":"back   ", bi, dist, dist>=rho-1e-4?"(EMPTY)":"(blocked)");
            if(surfside){ surf_sgn=sgn; surf_O[0]=O[0]; surf_O[1]=O[1]; surf_O[2]=O[2]; }
        }
        if(!surf_sgn){ printf("      VERDICT: no surface-consistent ball side (degenerate normal) -> skip\n"); blocked++; continue; }

        /* exact Delaunay arbiter on the surface side: is the nearest intruder inside
         * ABC's circumcircle? This is stable; the empty-ball metric is precision-
         * fragile at <1e-3 vox (dump is world float32, ULP ~5e-4 at coords ~4500). */
        int si=-1; double sbd=1e300;
        for(size_t i=0;i<nv;i++){ if((int)i==A||(int)i==B||(int)i==C) continue;
            double dx=V.p[i*3]-surf_O[0],dy=V.p[i*3+1]-surf_O[1],dz=V.p[i*3+2]-surf_O[2]; double dd=dx*dx+dy*dy+dz*dz;
            if(dd<sbd){sbd=dd;si=(int)i;} }
        double sdist=sqrt(sbd); int sempty=(sdist>=rho-1e-4);
        double pc[3]; v_sub(pc,&V.p[si*3],cc); double pheight=v_dot(pc,nu);
        double inplane=sqrt(fmax(0.0, v_dot(pc,pc)-pheight*pheight));
        double det=robust_incircle_value(a,b,c,&V.p[si*3],ns); int inside=(det>0.0);
        int she=0,shf=0;
        for(int t=0;t<3;t++) if(em_get(&em,si,trio[t])>0) she++;
        for(int t=0;t<3;t++){ int vtx=trio[t];
            for(size_t j=vf_off[vtx];j<vf_off[vtx+1];j++)
                for(size_t kk=vf_off[si];kk<vf_off[si+1];kk++)
                    if(vf_idx[j]==vf_idx[kk]) shf++; }
        printf("      surface intruder v%d: in-plane dist=%.4f vs rc=%.4f -> %s by %.4f vox (off-plane %.4f),"
               " shares %d/3 edges %d faces; exact in-circle=%+.2e -> %s\n",
               si, inplane, rc, inplane<rc?"INSIDE":"OUTSIDE", fabs(rc-inplane), pheight, she, shf,
               det, inside?"INSIDE circumcircle":"OUTSIDE circumcircle");

        /* verdict — the exact in-circle test is the arbiter */
        if(!simple){
            printf("      VERDICT: TRACER SKIP — boundary not a simple 3-cycle (bdeg!=2); fill_pinholes sets simple=0\n");
            tracer_skip++;
        } else if(inside){
            printf("      VERDICT: BLOCKED (cocircular) — adjacent fan apex v%d sits INSIDE ABC's circumcircle;\n"
                   "               filling ABC would overlap that fan / violate Delaunay. CORRECT skip.%s\n", si,
                   sempty ? "\n               (empty-ball read says empty, but by <1e-3 vox = dump rounding; the exact in-circle\n"
                            "               test is authoritative — internally the cube-local ball was non-empty.)" : "");
            blocked++;
        } else {
            printf("      VERDICT: FILLABLE (Delaunay-clear%s) yet it survived -> genuine internal miss:\n"
                   "               an edge stayed ES_FRONT (front never drained), so fill_pinholes never traced a\n"
                   "               closed 3-cycle here. Needs in-code fill_pinholes instrument to confirm.\n",
                   sempty?", empty ball":"");
            would_fill++;
        }
    }

    { long chi=(long)refV-(long)em.n+(long)nf;
      long genus=(2L*(long)comps-(long)nloops-chi)/2;
      printf("\n  topology: refV=%zu F=%zu E=%zu comps=%zu nonmanifold-E=%zu boundary-E=%zu loops=%d chi=%ld genus=%ld\n",
             refV, nf, em.n, comps, nmE, bE, nloops, chi, genus); }
    printf("  quality : mean-min-angle=%.1f deg  slivers(<10)=%ld (%.2f%%)  poor(<30)=%ld (%.1f%%)  mean-edge=%.3f vox\n"
           "            min-angle hist [<5 5-10 10-20 20-30 30-40 40-50 50-60]: %ld %ld %ld %ld %ld %ld %ld\n",
           minang_sum/(double)nf, nbad10, 100.0*(double)nbad10/(double)nf, nbad30, 100.0*(double)nbad30/(double)nf,
           elen_sum/(3.0*(double)nf), qhist[0],qhist[1],qhist[2],qhist[3],qhist[4],qhist[5],qhist[6]);
    printf("  summary : %d single-triangle holes | tracer-skip=%d  blocked/cocircular=%d  fillable-but-survived=%d\n",
           n3, tracer_skip, blocked, would_fill);
    free(V.p); free(Fc.p); free(vf_cnt); free(vf_off); free(vf_idx);
    free(em.key); free(em.val); free(bdeg); free(next); free(vis); free(uf);
    return 0;
}
