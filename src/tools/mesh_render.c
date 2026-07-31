/* mesh_render.c -- render one OBJ mesh to a high-resolution PNG for gallery QC.
 *
 * A self-contained software rasterizer (z-buffer, no GPU/GL). Adapted from
 * D:\repo\tools\objsequencetogif (camera/raster math) but: single high-res PNG
 * instead of a GIF sequence, HALF-LAMBERT lighting (wraps light around so the
 * whole surface reads, no black terminator), and surface colour driven by the
 * per-pixel NORMAL DIRECTION ("colours from descriptive angles") -- a clean
 * wrap shows a smooth colour gradient, a garbage/crumpled/fragmented sheet
 * shows chaotic speckle, which is exactly what we need to spot in a gallery.
 * Thin recto sheets are two-sided (no backface cull; the normal is flipped
 * toward the camera per pixel). View is a fixed 3/4 angle, camera auto-fit to
 * the mesh bbox, so every cube is framed the same way.
 *
 * Usage:
 *   mesh_render <in.obj> <out.png> [--size N] [--bg RRGGBB] [--gray] [--no-edges]
 *   mesh_render --selftest
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/* ---- tiny vector / matrix math ---- */
typedef struct { float x, y, z; } V3;
static V3 v3(float x, float y, float z){ V3 r; r.x=x; r.y=y; r.z=z; return r; }
static V3 vsub(V3 a, V3 b){ return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static V3 vadd(V3 a, V3 b){ return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static V3 vscl(V3 a, float s){ return v3(a.x*s, a.y*s, a.z*s); }
static float vdot(V3 a, V3 b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static V3 vcross(V3 a, V3 b){ return v3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static V3 vnorm(V3 a){ float l=sqrtf(vdot(a,a)); if(l<1e-8f) return v3(0,0,0); return vscl(a,1.0f/l); }
static float clamp01(float v){ return v<0.0f?0.0f:(v>1.0f?1.0f:v); }

typedef struct { float m[4][4]; } M4;
static M4 m4_id(void){ M4 r; memset(&r,0,sizeof r); r.m[0][0]=r.m[1][1]=r.m[2][2]=r.m[3][3]=1.0f; return r; }
static M4 m4_mul(M4 a, M4 b){
    M4 r; memset(&r,0,sizeof r);
    for(int i=0;i<4;i++) for(int j=0;j<4;j++){ float s=0; for(int k=0;k<4;k++) s+=a.m[i][k]*b.m[k][j]; r.m[i][j]=s; }
    return r;
}
/* transform point (w=1) with perspective divide */
static V3 xform_pt(M4 t, V3 p){
    float x=t.m[0][0]*p.x+t.m[0][1]*p.y+t.m[0][2]*p.z+t.m[0][3];
    float y=t.m[1][0]*p.x+t.m[1][1]*p.y+t.m[1][2]*p.z+t.m[1][3];
    float z=t.m[2][0]*p.x+t.m[2][1]*p.y+t.m[2][2]*p.z+t.m[2][3];
    float w=t.m[3][0]*p.x+t.m[3][1]*p.y+t.m[3][2]*p.z+t.m[3][3];
    if(fabsf(w)>1e-8f) return v3(x/w,y/w,z/w);
    return v3(x,y,z);
}

static M4 m4_lookat(V3 eye, V3 ctr, V3 up){
    V3 f=vnorm(vsub(ctr,eye)), r=vnorm(vcross(f,up)), u=vcross(r,f);
    M4 m=m4_id();
    m.m[0][0]=r.x; m.m[0][1]=r.y; m.m[0][2]=r.z;
    m.m[1][0]=u.x; m.m[1][1]=u.y; m.m[1][2]=u.z;
    m.m[2][0]=-f.x; m.m[2][1]=-f.y; m.m[2][2]=-f.z;
    m.m[0][3]=-vdot(r,eye); m.m[1][3]=-vdot(u,eye); m.m[2][3]=vdot(f,eye);
    return m;
}
static M4 m4_ortho(float l, float r, float b, float t, float n, float f){
    M4 m; memset(&m,0,sizeof m);
    m.m[0][0]=2.0f/(r-l); m.m[1][1]=2.0f/(t-b); m.m[2][2]=-2.0f/(f-n);
    m.m[0][3]=-(r+l)/(r-l); m.m[1][3]=-(t+b)/(t-b); m.m[2][3]=-(f+n)/(f-n); m.m[3][3]=1.0f;
    return m;
}
static M4 m4_viewport(int w, int h){
    M4 m; memset(&m,0,sizeof m);
    m.m[0][0]=w/2.0f; m.m[0][3]=w/2.0f;
    m.m[1][1]=-h/2.0f; m.m[1][3]=h/2.0f;   /* flip Y (screen down) */
    m.m[2][2]=0.5f; m.m[2][3]=0.5f; m.m[3][3]=1.0f;
    return m;
}

/* mvp + world camera-forward for an auto-fit ortho view from (azim, elev) degrees.
 * azim spins around the world Y axis, elev raises the camera above the XZ plane.
 * The up vector flips to world Z near the poles so a top-down view (elev ~ 90)
 * does not gimbal-lock against the default Y up. */
static void make_camera(V3 bmin, V3 bmax, int w, int h, float azim_deg,
                        float elev_deg, M4 *mvp, V3 *view_dir){
    V3 ctr=vscl(vadd(bmin,bmax),0.5f);
    V3 diag=vsub(bmax,bmin);
    float d=sqrtf(vdot(diag,diag)); if(d<1e-6f) d=1.0f;
    float az=azim_deg*3.14159265f/180.0f, el=elev_deg*3.14159265f/180.0f, dist=d*2.0f;
    V3 eye=vadd(ctr, v3(dist*cosf(el)*sinf(az), dist*sinf(el), dist*cosf(el)*cosf(az)));
    V3 up = (fabsf(elev_deg) > 80.0f) ? v3(0,0,1) : v3(0,1,0);
    M4 view=m4_lookat(eye,ctr,up);
    float radius=d*0.55f, aspect=(float)w/(float)h;
    M4 proj=m4_ortho(-radius*aspect, radius*aspect, -radius, radius, 0.01f, dist*4.0f);
    M4 vp=m4_viewport(w,h);
    *mvp=m4_mul(vp, m4_mul(proj, view));
    *view_dir=vnorm(v3(-view.m[2][0], -view.m[2][1], -view.m[2][2]));
}

/* ---- render one mesh into an RGBA buffer ---- */
static void render(const float *verts, size_t nv, const int32_t *faces, size_t nf,
                   const V3 *vn, const float *vcolors, int W, int H, int gray,
                   int edges, int ptsize, float azim, float elev,
                   unsigned char bg[3], unsigned char *img)
{
    /* bbox */
    V3 bmin=v3(1e30f,1e30f,1e30f), bmax=v3(-1e30f,-1e30f,-1e30f);
    for(size_t i=0;i<nv;i++){
        float x=verts[i*3+0],y=verts[i*3+1],z=verts[i*3+2];
        if(x<bmin.x)bmin.x=x; if(y<bmin.y)bmin.y=y; if(z<bmin.z)bmin.z=z;
        if(x>bmax.x)bmax.x=x; if(y>bmax.y)bmax.y=y; if(z>bmax.z)bmax.z=z;
    }
    M4 mvp; V3 vdir; make_camera(bmin,bmax,W,H,azim,elev,&mvp,&vdir);
    V3 light=vnorm(v3(0.5f,0.8f,0.6f));

    float *depth=(float*)malloc((size_t)W*H*sizeof(float));
    V3 *nbuf=(V3*)malloc((size_t)W*H*sizeof(V3));
    for(int i=0;i<W*H;i++){ depth[i]=1e30f; nbuf[i]=v3(0,0,0);
        img[i*4+0]=bg[0]; img[i*4+1]=bg[1]; img[i*4+2]=bg[2]; img[i*4+3]=255; }

    if(nf>0) for(size_t f=0; f<nf; f++){
        int32_t i0=faces[f*3+0], i1=faces[f*3+1], i2=faces[f*3+2];
        V3 v0=v3(verts[i0*3],verts[i0*3+1],verts[i0*3+2]);
        V3 v1=v3(verts[i1*3],verts[i1*3+1],verts[i1*3+2]);
        V3 v2=v3(verts[i2*3],verts[i2*3+1],verts[i2*3+2]);
        V3 s0=xform_pt(mvp,v0), s1=xform_pt(mvp,v1), s2=xform_pt(mvp,v2);
        int minX=(int)floorf(fminf(s0.x,fminf(s1.x,s2.x))); if(minX<0)minX=0;
        int maxX=(int)ceilf (fmaxf(s0.x,fmaxf(s1.x,s2.x))); if(maxX>W-1)maxX=W-1;
        int minY=(int)floorf(fminf(s0.y,fminf(s1.y,s2.y))); if(minY<0)minY=0;
        int maxY=(int)ceilf (fmaxf(s0.y,fmaxf(s1.y,s2.y))); if(maxY>H-1)maxY=H-1;
        /* edge functions for barycentric (2D) */
        float ax=s0.x,ay=s0.y,bx=s1.x,by=s1.y,cx=s2.x,cy=s2.y;
        float area=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
        if(fabsf(area)<1e-9f) continue;
        for(int y=minY;y<=maxY;y++) for(int x=minX;x<=maxX;x++){
            float px=x+0.5f, py=y+0.5f;
            float w0=((bx-px)*(cy-py)-(by-py)*(cx-px))/area;
            float w1=((cx-px)*(ay-py)-(cy-py)*(ax-px))/area;
            float w2=1.0f-w0-w1;
            if(w0<-1e-4f||w1<-1e-4f||w2<-1e-4f) continue;
            float z=w0*s0.z+w1*s1.z+w2*s2.z;
            int pi=y*W+x;
            if(z>=depth[pi]) continue;
            depth[pi]=z;
            V3 N=vnorm(vadd(vadd(vscl(vn[i0],w0),vscl(vn[i1],w1)),vscl(vn[i2],w2)));
            if(vdot(N,vdir)>0.0f) N=vscl(N,-1.0f);   /* two-sided: face camera */
            nbuf[pi]=N;
            float hl=vdot(N,light)*0.5f+0.5f; hl*=hl;  /* half-lambert (Valve) */
            /* base colour: per-vertex (heatmap) if given, else normal-direction */
            V3 base;
            if(vcolors) base=v3(
                vcolors[(size_t)i0*3+0]*w0+vcolors[(size_t)i1*3+0]*w1+vcolors[(size_t)i2*3+0]*w2,
                vcolors[(size_t)i0*3+1]*w0+vcolors[(size_t)i1*3+1]*w1+vcolors[(size_t)i2*3+1]*w2,
                vcolors[(size_t)i0*3+2]*w0+vcolors[(size_t)i1*3+2]*w1+vcolors[(size_t)i2*3+2]*w2);
            else base = gray ? v3(0.72f,0.72f,0.72f)
                             : v3(N.x*0.5f+0.5f, N.y*0.5f+0.5f, N.z*0.5f+0.5f);
            float amb=0.35f, k=amb+(1.0f-amb)*hl;
            img[pi*4+0]=(unsigned char)(clamp01(base.x*k)*255.0f);
            img[pi*4+1]=(unsigned char)(clamp01(base.y*k)*255.0f);
            img[pi*4+2]=(unsigned char)(clamp01(base.z*k)*255.0f);
        }
    }
    else if(nv>0){
        /* point cloud (no faces): splat each vertex as a small depth-tested disc,
         * then derive a screen-space normal from the depth buffer (clouds carry no
         * per-vertex normals) so a clean sheet reads smooth and a tangled/noisy
         * cloud reads as speckle -- same visual language as the mesh renders. The
         * depth ramp (near=warm, far=cool) cues 3D layout. */
        float dd=sqrtf(vdot(vsub(bmax,bmin),vsub(bmax,bmin))); if(dd<1e-6f) dd=1.0f;
        float wpp=(2.0f*dd*0.55f)/(float)H;       /* world units per pixel (ortho) */
        int R=ptsize<1?1:ptsize;
        float *vz=(float*)malloc((size_t)W*H*sizeof(float));
        for(int i=0;i<W*H;i++) vz[i]=0.0f;
        float vzmin=1e30f,vzmax=-1e30f;
        for(size_t i=0;i<nv;i++){
            V3 P=v3(verts[i*3],verts[i*3+1],verts[i*3+2]);
            V3 s=xform_pt(mvp,P);
            float pv=vdot(P,vdir);
            int cx=(int)lroundf(s.x), cy=(int)lroundf(s.y);
            for(int dy=-R;dy<=R;dy++) for(int dx=-R;dx<=R;dx++){
                if(dx*dx+dy*dy>R*R) continue;
                int x=cx+dx,y=cy+dy; if(x<0||x>=W||y<0||y>=H) continue;
                int pi=y*W+x;
                if(s.z<depth[pi]){ depth[pi]=s.z; vz[pi]=pv; }
            }
            if(pv<vzmin)vzmin=pv; if(pv>vzmax)vzmax=pv;
        }
        float vrange=vzmax-vzmin; if(vrange<1e-6f) vrange=1.0f;
        V3 plight=vnorm(v3(0.4f,0.5f,0.75f));
        for(int y=0;y<H;y++) for(int x=0;x<W;x++){
            int pi=y*W+x; if(depth[pi]>=1e29f) continue;
            float zc=vz[pi];
            float zl=(x>0   && depth[pi-1]<1e29f)?vz[pi-1]:zc;
            float zr=(x<W-1 && depth[pi+1]<1e29f)?vz[pi+1]:zc;
            float zu=(y>0   && depth[pi-W]<1e29f)?vz[pi-W]:zc;
            float zd=(y<H-1 && depth[pi+W]<1e29f)?vz[pi+W]:zc;
            V3 N=vnorm(v3(-(zr-zl)/(2.0f*wpp), -(zd-zu)/(2.0f*wpp), 1.0f));
            nbuf[pi]=N;
            float hl=vdot(N,plight)*0.5f+0.5f; hl*=hl;
            float t=(zc-vzmin)/vrange;
            V3 base = gray ? v3(0.72f,0.72f,0.72f)
                           : vadd(vscl(v3(0.95f,0.85f,0.55f),1.0f-t),
                                  vscl(v3(0.25f,0.40f,0.70f),t));
            float amb=0.35f,k=amb+(1.0f-amb)*hl;
            img[pi*4+0]=(unsigned char)(clamp01(base.x*k)*255.0f);
            img[pi*4+1]=(unsigned char)(clamp01(base.y*k)*255.0f);
            img[pi*4+2]=(unsigned char)(clamp01(base.z*k)*255.0f);
        }
        free(vz);
    }

    /* silhouette + crease edge overlay (dark), so sheet boundaries read crisply.
     * Meshes only: on a dense point splat, overlapping sheets make every pixel a
     * depth discontinuity, so the overlay would paint the whole frame dark. */
    if(edges && nf>0){
        for(int y=1;y<H-1;y++) for(int x=1;x<W-1;x++){
            int idx=y*W+x; if(depth[idx]>=1e29f) continue;
            int is_edge=0; float dc=depth[idx]; V3 nc=nbuf[idx];
            for(int dy=-1;dy<=1&&!is_edge;dy++) for(int dx=-1;dx<=1&&!is_edge;dx++){
                if(!dx&&!dy) continue; int n=(y+dy)*W+(x+dx);
                if(depth[n]>=1e29f){ is_edge=1; }
                else if(fabsf(depth[n]-dc)>0.01f){ is_edge=1; }
                else if(vdot(nc,nbuf[n])<0.55f){ is_edge=1; }
            }
            if(is_edge){ img[idx*4+0]=30; img[idx*4+1]=30; img[idx*4+2]=30; }
        }
    }
    free(depth); free(nbuf);
}

/* area-weighted per-vertex normals */
static void compute_normals(const float *verts, size_t nv,
                            const int32_t *faces, size_t nf, V3 *vn){
    for(size_t i=0;i<nv;i++) vn[i]=v3(0,0,0);
    for(size_t f=0;f<nf;f++){
        int32_t a=faces[f*3],b=faces[f*3+1],c=faces[f*3+2];
        V3 va=v3(verts[a*3],verts[a*3+1],verts[a*3+2]);
        V3 vb=v3(verts[b*3],verts[b*3+1],verts[b*3+2]);
        V3 vc=v3(verts[c*3],verts[c*3+1],verts[c*3+2]);
        V3 fn=vcross(vsub(vb,va),vsub(vc,va));   /* area-weighted (unnormalized) */
        vn[a]=vadd(vn[a],fn); vn[b]=vadd(vn[b],fn); vn[c]=vadd(vn[c],fn);
    }
    for(size_t i=0;i<nv;i++) vn[i]=vnorm(vn[i]);
}

static int run_selftest(void){
    /* a single triangle must paint some non-background pixels */
    float verts[9]={0,0,0, 1,0,0, 0,1,0};
    int32_t faces[3]={0,1,2};
    int W=64,H=64; V3 vn[3];
    compute_normals(verts,3,faces,1,vn);
    unsigned char bg[3]={255,255,255};
    unsigned char *img=(unsigned char*)malloc((size_t)W*H*4);
    render(verts,3,faces,1,vn,NULL,W,H,0,1,2,35.0f,25.0f,bg,img);
    long painted=0; for(int i=0;i<W*H;i++) if(img[i*4]!=255||img[i*4+1]!=255||img[i*4+2]!=255) painted++;
    /* a point cloud (no faces) must also paint pixels via the splat path */
    float pts[12]={0,0,0, 1,0,0, 0,1,0, 1,1,0};
    render(pts,4,NULL,0,NULL,NULL,W,H,0,0,3,35.0f,25.0f,bg,img);
    long ppx=0; for(int i=0;i<W*H;i++) if(img[i*4]!=255||img[i*4+1]!=255||img[i*4+2]!=255) ppx++;
    free(img);
    int fail = (painted < 50) || (ppx < 4);
    fprintf(stderr,"[selftest] triangle painted %ld px, points painted %ld px -> %s\n",
            painted, ppx, fail?"FAIL":"ok");
    fprintf(stderr,"=== mesh_render selftest %s ===\n", fail?"FAILED":"PASSED");
    return fail?1:0;
}

static void parse_hex(const char *s, unsigned char out[3]){
    unsigned v=(unsigned)strtoul(s,NULL,16);
    out[0]=(v>>16)&0xFF; out[1]=(v>>8)&0xFF; out[2]=v&0xFF;
}

/* Parse per-vertex RGB from an extended OBJ ("v x y z r g b" -- the format
 * ObjIO_write_per_vertex_color emits, which ObjIO_read silently drops). Returns
 * a malloc'd [nv*3] array (caller frees) or NULL if the file carries no vertex
 * colours. Parses v-lines in file order, matching ObjIO_read's vertex order. */
static float *parse_vertex_colors(const char *path, size_t nv){
    FILE *fp=fopen(path,"rb"); if(!fp) return NULL;
    float *col=(float*)malloc((nv?nv:1)*3*sizeof(float));
    size_t vi=0; int any=0; char line[512];
    while(fgets(line,sizeof line,fp)){
        if(line[0]!='v'||line[1]!=' ') continue;
        float x,y,z,r,g,b;
        int n=sscanf(line+2,"%f %f %f %f %f %f",&x,&y,&z,&r,&g,&b);
        if(vi<nv){
            if(n>=6){ col[vi*3+0]=r; col[vi*3+1]=g; col[vi*3+2]=b; any=1; }
            else    { col[vi*3+0]=col[vi*3+1]=col[vi*3+2]=0.72f; }
        }
        vi++;
    }
    fclose(fp);
    if(!any){ free(col); return NULL; }
    return col;
}

/* Blit a WxH RGBA tile into a (cols*W)x(rows*H) composite at grid cell (gx,gy). */
static void blit_tile(const unsigned char *tile, int W, int H,
                      unsigned char *comp, int cols, int gx, int gy){
    int CW=cols*W;
    for(int y=0;y<H;y++) for(int x=0;x<W;x++){
        int s=(y*W+x)*4, d=(((gy*H)+y)*CW+(gx*W)+x)*4;
        comp[d+0]=tile[s+0]; comp[d+1]=tile[s+1]; comp[d+2]=tile[s+2]; comp[d+3]=tile[s+3];
    }
}

int main(int argc, char **argv){
    if(argc>=2 && !strcmp(argv[1],"--selftest")) return run_selftest();
    if(argc<3){
        fprintf(stderr,"Usage: %s <in.obj> <out.png> [--size N] [--bg RRGGBB] [--gray] [--no-edges] [--ptsize N]\n"
                       "                              [--azim DEG] [--elev DEG] [--views] [--vcolor]\n"
                       "       %s --selftest\n"
                       "  --azim/--elev : camera angle (default 35/25 deg)\n"
                       "  --views       : 2x2 contact sheet (TL iso, TR top, BL front, BR side); --size = per-tile\n"
                       "  --vcolor      : shade per-vertex OBJ colours (e.g. a developability heatmap) not normals\n"
                       "  (a faceless OBJ is rendered as a splatted point cloud; --ptsize sets splat radius)\n",
                       argv[0], argv[0]);
        return 1;
    }
    const char *in=argv[1], *out=argv[2];
    int W=768,H=768,gray=0,edges=1,ptsize=2,views=0,vcolor=0;
    float azim=35.0f, elev=25.0f;
    unsigned char bg[3]={255,255,255};
    for(int i=3;i<argc;i++){
        if(!strcmp(argv[i],"--size")&&i+1<argc){ W=H=atoi(argv[++i]); }
        else if(!strcmp(argv[i],"--bg")&&i+1<argc){ parse_hex(argv[++i],bg); }
        else if(!strcmp(argv[i],"--gray")){ gray=1; }
        else if(!strcmp(argv[i],"--no-edges")){ edges=0; }
        else if(!strcmp(argv[i],"--ptsize")&&i+1<argc){ ptsize=atoi(argv[++i]); }
        else if(!strcmp(argv[i],"--azim")&&i+1<argc){ azim=(float)atof(argv[++i]); }
        else if(!strcmp(argv[i],"--elev")&&i+1<argc){ elev=(float)atof(argv[++i]); }
        else if(!strcmp(argv[i],"--views")){ views=1; }
        else if(!strcmp(argv[i],"--vcolor")||!strcmp(argv[i],"--vcolour")){ vcolor=1; }
        else { fprintf(stderr,"unknown arg %s\n",argv[i]); return 1; }
    }
    if(W<16||W>8192){ fprintf(stderr,"--size out of range\n"); return 1; }
    if(ptsize<1)ptsize=1; if(ptsize>16)ptsize=16;

    Arena_T arena=Arena_new();
    float *verts=NULL; int32_t *faces=NULL; size_t nv=0,nf=0;
    if(ObjIO_read(arena,in,&verts,&nv,&faces,&nf)!=0){
        fprintf(stderr,"mesh_render: cannot read %s\n",in); Arena_dispose(&arena); return 1;
    }
    float *vcolors=NULL;
    if(vcolor){
        vcolors=parse_vertex_colors(in,nv);
        if(!vcolors) fprintf(stderr,"  (--vcolor: %s has no per-vertex colours; using normals)\n", in);
    }
    V3 *vn=NULL;
    if(nv>0){ vn=(V3*)ARENA_ALLOC(arena,(long)(nv*sizeof(V3)));
              if(nf>0) compute_normals(verts,nv,faces,nf,vn); }

    /* --views: 2x2 contact sheet -- iso / top / front / side -- so an edge-on
     * angle reveals fold-vs-bridge that the single 3/4 view can't. Otherwise one
     * view from (azim,elev). --size is the PER-TILE size in --views mode. */
    int ok=0;
    ves_ensure_parent_dir(out);
    if(views){
        const float vw[4][2]={{35,25},{0,89},{0,0},{90,0}};  /* iso, top, front, side */
        int CW=2*W, CH=2*H;
        unsigned char *comp=(unsigned char*)malloc((size_t)CW*CH*4);
        unsigned char *tile=(unsigned char*)malloc((size_t)W*H*4);
        if(!comp||!tile){ fprintf(stderr,"oom\n"); free(comp); free(tile); Arena_dispose(&arena); return 1; }
        for(int k=0;k<4;k++){
            if(nv==0) for(int i=0;i<W*H;i++){ tile[i*4]=bg[0];tile[i*4+1]=bg[1];tile[i*4+2]=bg[2];tile[i*4+3]=255; }
            else render(verts,nv,faces,nf,vn,vcolors,W,H,gray,edges,ptsize,vw[k][0],vw[k][1],bg,tile);
            blit_tile(tile,W,H,comp,2,k%2,k/2);
        }
        ok=stbi_write_png(out,CW,CH,4,comp,CW*4);
        free(comp); free(tile);
        fprintf(stderr,"mesh_render: %zu v %zu f -> %dx%d %s [views: TL iso, TR top, BL front, BR side]\n",
                nv,nf,CW,CH,out);
    } else {
        unsigned char *img=(unsigned char*)malloc((size_t)W*H*4);
        if(!img){ fprintf(stderr,"oom\n"); Arena_dispose(&arena); free(vcolors); return 1; }
        if(nv==0) for(int i=0;i<W*H;i++){ img[i*4]=bg[0];img[i*4+1]=bg[1];img[i*4+2]=bg[2];img[i*4+3]=255; }
        else render(verts,nv,faces,nf,vn,vcolors,W,H,gray,edges,ptsize,azim,elev,bg,img);
        ok=stbi_write_png(out,W,H,4,img,W*4);
        free(img);
        fprintf(stderr,"mesh_render: %zu v %zu f -> %dx%d %s (azim=%.0f elev=%.0f%s)\n",
                nv,nf,W,H,out,(double)azim,(double)elev,vcolors?" vcolor":"");
    }
    free(vcolors);
    Arena_dispose(&arena);
    if(!ok){ fprintf(stderr,"mesh_render: PNG write failed %s\n",out); return 1; }
    return 0;
}
