/* dev_metric.c -- discrete angle-defect Gaussian curvature (developability).
 * See dev_metric.h. Cold paths only (malloc scratch). */
#include "dev_metric.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEVM_PI 3.14159265358979323846

static double clampd(double v, double lo, double hi){ return v<lo?lo:(v>hi?hi:v); }

static uint64_t ekey(int32_t a, int32_t b){
    uint32_t lo = (uint32_t)(a<b?a:b), hi = (uint32_t)(a<b?b:a);
    return ((uint64_t)lo<<32) | hi;
}
static int cmp_u64(const void *pa, const void *pb){
    uint64_t a=*(const uint64_t*)pa, b=*(const uint64_t*)pb;
    return (a<b)?-1:((a>b)?1:0);
}

/* corner angle at p0 in triangle (p0,p1,p2); 0 if degenerate */
static double corner_angle(const float *p0, const float *p1, const float *p2){
    double u0=p1[0]-p0[0], u1=p1[1]-p0[1], u2=p1[2]-p0[2];
    double v0=p2[0]-p0[0], v1=p2[1]-p0[1], v2=p2[2]-p0[2];
    double lu=sqrt(u0*u0+u1*u1+u2*u2), lv=sqrt(v0*v0+v1*v1+v2*v2);
    if(lu<1e-12 || lv<1e-12) return 0.0;
    double c=(u0*v0+u1*v1+u2*v2)/(lu*lv);
    return acos(clampd(c,-1.0,1.0));
}

void DevMetric_compute(const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       double *defect, unsigned char *boundary,
                       unsigned char *used)
{
    double *angsum = (double*)calloc(nv?nv:1, sizeof(double));
    memset(boundary, 0, nv);
    memset(used, 0, nv);

    /* boundary: an undirected edge used exactly once is a boundary edge. */
    size_t ne = 3*nf;
    uint64_t *keys = (uint64_t*)malloc((ne?ne:1)*sizeof(uint64_t));
    size_t k=0;
    for(size_t f=0; f<nf; f++){
        int32_t a=faces[f*3+0], b=faces[f*3+1], c=faces[f*3+2];
        keys[k++]=ekey(a,b); keys[k++]=ekey(b,c); keys[k++]=ekey(c,a);
    }
    qsort(keys, ne, sizeof(uint64_t), cmp_u64);
    for(size_t i=0; i<ne; ){
        size_t j=i+1; while(j<ne && keys[j]==keys[i]) j++;
        if(j-i==1){
            int32_t lo=(int32_t)(keys[i]>>32), hi=(int32_t)(keys[i]&0xffffffffu);
            boundary[lo]=1; boundary[hi]=1;
        }
        i=j;
    }
    free(keys);

    for(size_t f=0; f<nf; f++){
        int32_t i0=faces[f*3+0], i1=faces[f*3+1], i2=faces[f*3+2];
        const float *p0=&verts[(size_t)i0*3], *p1=&verts[(size_t)i1*3], *p2=&verts[(size_t)i2*3];
        angsum[i0]+=corner_angle(p0,p1,p2);
        angsum[i1]+=corner_angle(p1,p2,p0);
        angsum[i2]+=corner_angle(p2,p0,p1);
        used[i0]=used[i1]=used[i2]=1;
    }
    for(size_t i=0;i<nv;i++){
        if(!used[i]){ defect[i]=0.0; continue; }
        double datum = boundary[i] ? DEVM_PI : 2.0*DEVM_PI;
        defect[i] = datum - angsum[i];
    }
    free(angsum);
}

void DevMetric_interior_stats(const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              double k_thresh, DevMetricStats *out)
{
    memset(out, 0, sizeof(*out));
    if(nv==0 || nf==0) return;
    double *defect=(double*)malloc(nv*sizeof(double));
    unsigned char *boundary=(unsigned char*)malloc(nv);
    unsigned char *used=(unsigned char*)malloc(nv);
    DevMetric_compute(verts,nv,faces,nf,defect,boundary,used);
    for(size_t i=0;i<nv;i++){
        if(!used[i]) continue;
        if(boundary[i]){ out->n_boundary++; continue; }
        double ad=fabs(defect[i]);
        out->n_interior++;
        out->sum_abs_k_interior += ad;
        if(ad>out->max_abs_k) out->max_abs_k=ad;
        if(ad>k_thresh) out->n_interior_bad++;
    }
    out->frac_interior_bad = out->n_interior ?
        (double)out->n_interior_bad/(double)out->n_interior : 0.0;
    free(defect); free(boundary); free(used);
}
