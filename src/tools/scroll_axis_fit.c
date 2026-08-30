/* scroll_axis_fit.c -- estimate the scroll central axis (umbilicus line) and the
 * radial wrap-spacing from one or more welded-region OBJs.
 *
 * A 640^3 sub-block can't recover the scroll's winding frame (its local spiral
 * fit collapses), but a larger region that brackets the umbilicus can: PCA over
 * a roughly-concentric wrap surface gives the axis DIRECTION (the "odd one out"
 * eigenvector, PCA_scroll_axis) and a CENTROID ~ the umbilicus. The radial
 * distribution of vertices about that line is multi-modal with one peak per
 * wrap, so the median peak spacing estimates the radial gain (vox per turn).
 *
 * Prints the axis line + wrap-spacing to feed:
 *   scroll_atlas ... --axis-point z y x --axis-dir z y x --wrap-spacing B
 *
 * Usage:
 *   scroll_axis_fit <region1.obj> [region2.obj ...] [--binw W] [--umbilicus-yx Y X]
 *   scroll_axis_fit --selftest
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"
#include "../common/pca.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* perpendicular distance of point p (z,y,x) from the line through cen along
 * unit direction a. */
static double radius_from_axis(const float *p, const double cen[3], const double a[3])
{
    double d0 = (double)p[0]-cen[0], d1 = (double)p[1]-cen[1], d2 = (double)p[2]-cen[2];
    double t  = d0*a[0] + d1*a[1] + d2*a[2];
    double q0 = d0 - t*a[0], q1 = d1 - t*a[1], q2 = d2 - t*a[2];
    return sqrt(q0*q0 + q1*q1 + q2*q2);
}

static int cmp_double(const void *a, const void *b)
{ double x=*(const double*)a, y=*(const double*)b; return (x<y)?-1:(x>y)?1:0; }

/* radii[n] (any order). Histogram (bin width binw), light smooth, find local-max
 * peaks above 25% of the max bin, return the median spacing between consecutive
 * peak radii. Reports peaks if verbose. Returns <=0 if <2 peaks. */
static double estimate_wrap_spacing(const double *radii, size_t n, double binw,
                                    int verbose)
{
    if (n < 2 || binw <= 0) return -1.0;
    double rmin=1e300, rmax=-1e300;
    for (size_t i=0;i<n;i++){ if(radii[i]<rmin)rmin=radii[i]; if(radii[i]>rmax)rmax=radii[i]; }
    if (rmax<=rmin) return -1.0;
    size_t nb = (size_t)((rmax-rmin)/binw) + 1;
    if (nb < 3) return -1.0;
    double *h = (double*)calloc(nb, sizeof(double));
    for (size_t i=0;i<n;i++){ size_t b=(size_t)((radii[i]-rmin)/binw); if(b>=nb)b=nb-1; h[b]+=1.0; }
    /* moving-average smooth (half-window 1) */
    double *s = (double*)malloc(nb*sizeof(double));
    for (size_t b=0;b<nb;b++){ double a=0; int c=0; for(long w=(long)b-1;w<=(long)b+1;w++){ if(w<0||w>=(long)nb)continue; a+=h[w]; c++; } s[b]=a/c; }
    double hmax=0; for(size_t b=0;b<nb;b++) if(s[b]>hmax) hmax=s[b];
    double thr = 0.25*hmax;
    /* peak radii: local maxima above threshold (>= both neighbors) */
    double *peaks = (double*)malloc(nb*sizeof(double)); size_t np=0;
    for (size_t b=0;b<nb;b++){
        double l = (b>0)?s[b-1]:0.0, rr = (b+1<nb)?s[b+1]:0.0;
        if (s[b]>=thr && s[b]>=l && s[b]>=rr) {
            double rc = rmin + ((double)b+0.5)*binw;
            /* merge plateau peaks closer than binw */
            if (np>0 && rc-peaks[np-1] < binw*1.5) peaks[np-1]=0.5*(peaks[np-1]+rc);
            else peaks[np++]=rc;
        }
    }
    double b_est = -1.0;
    if (np>=2){
        double *d=(double*)malloc((np-1)*sizeof(double));
        for(size_t i=0;i+1<np;i++) d[i]=peaks[i+1]-peaks[i];
        qsort(d,np-1,sizeof(double),cmp_double);
        b_est = d[(np-1)/2];
        free(d);
    }
    if (verbose){
        fprintf(stderr,"  radii: n=%zu range=[%.1f,%.1f] bins=%zu peaks=%zu\n", n,rmin,rmax,nb,np);
        if (np){ fprintf(stderr,"  peak radii:"); for(size_t i=0;i<np && i<40;i++) fprintf(stderr," %.1f",peaks[i]); fprintf(stderr,"\n"); }
    }
    free(h); free(s); free(peaks);
    return b_est;
}

/* concat all verts from the OBJs (faces ignored). malloc'd; caller frees. */
static float *load_clouds(char **paths, int n, size_t *out_nv)
{
    float *V=NULL; size_t nv=0, cap=0;
    for (int a=0;a<n;a++){
        Arena_T ar=Arena_new();
        float *v=NULL; int32_t *f=NULL; size_t lnv=0,lnf=0;
        if (ObjIO_read(ar,paths[a],&v,&lnv,&f,&lnf)==0 && lnv){
            if (nv+lnv>cap){ cap=(nv+lnv)*2+1024; V=(float*)realloc(V,cap*3*sizeof(float)); }
            memcpy(V+nv*3, v, lnv*3*sizeof(float)); nv+=lnv;
            fprintf(stderr,"  loaded %s: %zu v\n", paths[a], lnv);
        } else fprintf(stderr,"  skip (unreadable/empty) %s\n", paths[a]);
        Arena_dispose(&ar);
    }
    *out_nv=nv; return V;
}

static void fit_and_report(const float *V, size_t nv, double binw)
{
    float axf[3]={0,0,0}, cenf[3]={0,0,0};
    PCA_scroll_axis(V, nv, axf, cenf);
    double a[3]={axf[0],axf[1],axf[2]}, cen[3]={cenf[0],cenf[1],cenf[2]};
    double nrm=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); if(nrm<1e-12)nrm=1.0;
    a[0]/=nrm; a[1]/=nrm; a[2]/=nrm;
    double *rad=(double*)malloc(nv*sizeof(double));
    for (size_t i=0;i<nv;i++) rad[i]=radius_from_axis(&V[i*3], cen, a);
    double b=estimate_wrap_spacing(rad, nv, binw, 1);
    free(rad);
    printf("axis_point (z y x): %.2f %.2f %.2f\n", cen[0],cen[1],cen[2]);
    printf("axis_dir   (z y x): %.4f %.4f %.4f\n", a[0],a[1],a[2]);
    printf("wrap_spacing (vox/turn): %.3f\n", b);
    printf("\nscroll_atlas flags:\n  --axis-point %.2f %.2f %.2f --axis-dir %.4f %.4f %.4f --wrap-spacing %.3f\n",
           cen[0],cen[1],cen[2], a[0],a[1],a[2], b>0?b:1.0);
}

static int run_selftest(void)
{
    /* concentric cylinders around the Z axis (axis=(1,0,0) in z,y,x) centered at
     * (y,x)=(40,60); radii 5,10,15,20 (spacing 5); modest height in z. */
    int fail=0;
    double cy=40, cx=60; int RAD[4]={5,10,15,20}, NT=120, NZ=12;
    size_t nv=(size_t)4*NT*NZ;
    float *V=(float*)malloc(nv*3*sizeof(float)); size_t k=0;
    for (int ri=0;ri<4;ri++) for(int ti=0;ti<NT;ti++) for(int zi=0;zi<NZ;zi++){
        double th=2.0*M_PI*ti/NT;
        V[k*3+0]=(float)zi;                       /* z (axial) */
        V[k*3+1]=(float)(cy+RAD[ri]*cos(th));     /* y */
        V[k*3+2]=(float)(cx+RAD[ri]*sin(th));     /* x */
        k++;
    }
    float axf[3], cenf[3]; PCA_scroll_axis(V,nv,axf,cenf);
    double az=fabs(axf[0]);   /* axis should be ~ +/- Z */
    int f1=(az<0.9); fail|=f1;
    fprintf(stderr,"[selftest] axis=(%.3f,%.3f,%.3f) |z|=%.3f want~1 -> %s\n",axf[0],axf[1],axf[2],az,f1?"FAIL":"ok");
    double a[3]={axf[0],axf[1],axf[2]}; double cen[3]={cenf[0],cenf[1],cenf[2]};
    double *rad=(double*)malloc(nv*sizeof(double));
    for(size_t i=0;i<nv;i++) rad[i]=radius_from_axis(&V[i*3],cen,a);
    double b=estimate_wrap_spacing(rad,nv,1.0,0);
    int f2=!(b>4.0 && b<6.0); fail|=f2;   /* expect ~5 */
    fprintf(stderr,"[selftest] wrap_spacing=%.3f want~5 -> %s\n", b, f2?"FAIL":"ok");
    /* centroid should sit near the umbilicus (y=40,x=60) */
    int f3=!(fabs(cenf[1]-cy)<2.0 && fabs(cenf[2]-cx)<2.0); fail|=f3;
    fprintf(stderr,"[selftest] centroid (y,x)=(%.1f,%.1f) want(40,60) -> %s\n",cenf[1],cenf[2],f3?"FAIL":"ok");
    free(rad); free(V);
    fprintf(stderr,"=== scroll_axis_fit selftest %s ===\n", fail?"FAILED":"PASSED");
    return fail?1:0;
}

int main(int argc, char **argv)
{
    if (argc>=2 && !strcmp(argv[1],"--selftest")) return run_selftest();
    if (argc<2){
        fprintf(stderr,"Usage: %s <region1.obj> [region2.obj ...] [--binw W]\n"
                       "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }
    double binw=1.0;
    char *paths[256]; int np=0;
    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"--binw") && i+1<argc) binw=atof(argv[++i]);
        else if (np<256) paths[np++]=argv[i];
    }
    if (np<1){ fprintf(stderr,"no input OBJs\n"); return 1; }
    size_t nv=0; float *V=load_clouds(paths, np, &nv);
    if (nv<16){ fprintf(stderr,"too few verts (%zu)\n", nv); free(V); return 1; }
    fprintf(stderr,"scroll_axis_fit: %d region(s), %zu verts\n", np, nv);
    fit_and_report(V, nv, binw);
    free(V);
    return 0;
}
