#include "asm_metric.h"
#include "../common/uv_guard.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int AsmMetric_prepare(const AsmChart *c, AsmMetricFace *ref)
{
    if (!c || !c->nf || !c->nv || !c->xyz || !c->faces || !ref) return -1;
    for (size_t i = 0; i < c->nf; i++) {
        const int32_t *f = c->faces+3*i;
        for (size_t k = 0; k < 3; k++) if (f[k] < 0 || (size_t)f[k] >= c->nv) return -1;
        double e[3], d[3], cross[3], length2 = 0, dot = 0;
        for (size_t k = 0; k < 3; k++) {
            e[k] = (double)c->xyz[3*(size_t)f[1]+k]-c->xyz[3*(size_t)f[0]+k];
            d[k] = (double)c->xyz[3*(size_t)f[2]+k]-c->xyz[3*(size_t)f[0]+k];
            length2 += e[k]*e[k]; dot += e[k]*d[k];
        }
        for (size_t k = 0; k < 3; k++) cross[k] = e[(k+1)%3]*d[(k+2)%3]-e[(k+2)%3]*d[(k+1)%3];
        double length = sqrt(length2), area2 = hypot(hypot(cross[0],cross[1]),cross[2]);
        if (!(length > 0) || !(area2 > 0) || !isfinite(length) || !isfinite(area2) || !isfinite(dot)) return -1;
        ref[i] = (AsmMetricFace){length,dot/length,area2/length,.5*area2};
    }
    return 0;
}

static double am_uv(const void *uv, size_t i, int precise)
{
    return precise ? ((const double *)uv)[i] : (double)((const float *)uv)[i];
}

static void am_measure(const AsmChart *c, const AsmMetricFace *ref,
                        const void *uv, int precise, int orientation, AsmMetricStats *out)
{
    memset(out,0,sizeof *out); out->minimum = DBL_MAX;
    for (size_t i = 0; i < c->nf; i++) {
        const int32_t *f = c->faces+3*i; const AsmMetricFace *r = ref+i;
        double u0 = am_uv(uv,2*(size_t)f[0],precise), v0 = am_uv(uv,2*(size_t)f[0]+1,precise);
        double a = (am_uv(uv,2*(size_t)f[1],precise)-u0)/r->length;
        double c0 = (am_uv(uv,2*(size_t)f[1]+1,precise)-v0)/r->length;
        double b = (am_uv(uv,2*(size_t)f[2],precise)-u0-r->along*a)/r->height;
        double d = (am_uv(uv,2*(size_t)f[2]+1,precise)-v0-r->along*c0)/r->height;
        double det = a*d-b*c0;
        double hi = .5*(hypot(a+d,c0-b)+hypot(a-d,c0+b));
        double lo = hi > 0 ? fabs(det)/hi : 0;
        out->area += r->area;
        if (!(orientation*det > 0) || !(lo > 0) || !isfinite(lo) || !isfinite(hi)) { out->invalid++; continue; }
        out->minimum = fmin(out->minimum,lo); out->maximum = fmax(out->maximum,hi);
        if (AsmMetric_within10(lo,hi)) out->within10 += r->area;
        if (AsmMetric_within25(lo,hi)) out->within25 += r->area;
    }
}

void AsmMetric_measure(const AsmChart *c, const AsmMetricFace *ref,
                       const float *uv, AsmMetricStats *out)
{
    am_measure(c,ref,uv,0,1,out);
}

void AsmMetric_measure_encoded(const AsmChart *c, const AsmMetricFace *ref,
                               const float *uv, AsmMetricStats *out)
{
    am_measure(c,ref,uv,0,(c->flags & ASM_CHART_MIRROR) ? -1 : 1,out);
}

void AsmMetric_measure_encoded_uv64(const AsmChart *c, const AsmMetricFace *ref,
                                    const double *uv, AsmMetricStats *out)
{
    am_measure(c,ref,uv,1,(c->flags & ASM_CHART_MIRROR) ? -1 : 1,out);
}

int AsmMetric_preserved(const AsmMetricStats *before, const AsmMetricStats *after)
{
    if (!(before->area > 0) || !isfinite(before->area) || before->area != after->area ||
        before->invalid || after->invalid) return 0;
    /* This only absorbs summation roundoff; it cannot hide an original face. */
    double tolerance = 32*DBL_EPSILON*before->area;
    return after->within10+tolerance >= fmin(.95*before->area,before->within10) &&
           after->within25+tolerance >= fmin(.99*before->area,before->within25);
}

int AsmMetric_local_face_preserved(const AsmChart *c, size_t t,
                                  const AsmMetricFace *ref, const float *trial)
{
    const int32_t *f=c->faces+3*t;int32_t packed_face[]={0,1,2};
    double old[6],next[6];
    for(int k=0;k<3;k++)for(int d=0;d<2;d++){
        old[2*k+d]=c->uv[2*(size_t)f[k]+d];next[2*k+d]=trial[2*(size_t)f[k]+d];
    }
    if(!UvGuard_interval(old,next,packed_face,1) || !UvGuard_positive(next,packed_face))return 0;
    AsmChart triangle=*c;triangle.faces=f;triangle.nf=1;
    AsmMetricStats before,after;AsmMetric_measure(&triangle,ref,c->uv,&before);AsmMetric_measure(&triangle,ref,trial,&after);
    if(before.invalid || after.invalid)return 0;
    return !AsmMetric_within25(before.minimum,before.maximum) || AsmMetric_within25(after.minimum,after.maximum);
}

int AsmMetric_encoded_faces_preserved(const AsmChart *c, const AsmMetricFace *ref, const double *encoded)
{
    return AsmMetric_encoded_faces_preserved_parity(c, ref, encoded, (c->flags & ASM_CHART_MIRROR) != 0);
}

int AsmMetric_encoded_faces_preserved_parity(const AsmChart *c, const AsmMetricFace *ref, const double *encoded,
                                             int encoded_mirror)
{
    for (size_t t = 0; t < c->nf; t++) {
        double lo[2], hi[2];
        const int32_t *f = c->faces+3*t; const AsmMetricFace *r = ref+t;
        for (int side = 0; side < 2; side++) {
            const void *uv = side ? (const void *)encoded : c->placed_uv ? (const void *)c->placed_uv : (const void *)c->uv;
            int precise = side || c->placed_uv != NULL;
            double u = am_uv(uv,2*(size_t)f[0],precise), v = am_uv(uv,2*(size_t)f[0]+1,precise);
            double a = (am_uv(uv,2*(size_t)f[1],precise)-u)/r->length;
            double b = (am_uv(uv,2*(size_t)f[1]+1,precise)-v)/r->length;
            double x = (am_uv(uv,2*(size_t)f[2],precise)-u-r->along*a)/r->height;
            double y = (am_uv(uv,2*(size_t)f[2]+1,precise)-v-r->along*b)/r->height;
            double det = a*y-b*x;
            int orientation = side ? (encoded_mirror ? -1 : 1) : (precise && (c->flags & ASM_CHART_MIRROR) ? -1 : 1);
            hi[side] = .5*(hypot(a+y,b-x)+hypot(a-y,b+x)); lo[side] = hi[side] > 0 ? fabs(det)/hi[side] : 0;
            if (!(orientation*det > 0) || !(lo[side] > 0) || !isfinite(hi[side])) return 0;
        }
        if (AsmMetric_within25(lo[0],hi[0]) && !AsmMetric_within25(lo[1],hi[1])) return 0;
    } return 1;
}

int AsmMetric_selftest(void)
{
    float xyz[] = {0,0,0, 10,0,0, 0,10,0}, uv[] = {0,0, 10,0, 0,10};
    int32_t face[] = {0,1,2}; AsmChart c = {0}; AsmMetricFace ref;
    c.nv = 3; c.nf = 1; c.xyz = xyz; c.uv = uv; c.faces = face;
    AsmMetricStats original, trial; int fails = 0;
    if (!AsmMetric_within10(1+DBL_EPSILON,1) ||
        !AsmMetric_within10(.9-2e-11,1.1+2e-11) ||
        !AsmMetric_within25(.75-2e-11,1.25+2e-11) ||
        AsmMetric_within10(.9-2*ASM_METRIC_ROUNDOFF,1) ||
        AsmMetric_within10(1,1.1+2*ASM_METRIC_ROUNDOFF) ||
        AsmMetric_within25(.75-2*ASM_METRIC_ROUNDOFF,1) ||
        AsmMetric_within25(1,1.25+2*ASM_METRIC_ROUNDOFF) ||
        AsmMetric_within25(NAN,1) || AsmMetric_within25(1,INFINITY)) fails++;
    if (AsmMetric_prepare(&c,&ref)) fails++;
    AsmMetric_measure(&c,&ref,uv,&original);
    if (original.area != 50 || original.within10 != 50 || original.minimum != 1 || original.maximum != 1) fails++;
    uv[2] = 8; AsmMetric_measure(&c,&ref,uv,&trial);
    if (AsmMetric_preserved(&original,&trial) || trial.within10 || trial.within25 != 50) fails++;
    if (!AsmMetric_preserved(&trial,&original) || !AsmMetric_preserved(&trial,&trial)) fails++;
    /* A translation and rotation preserve metric even in an absolute frame. */
    float rotated[] = {30000,400, 30000,410, 29990,400};
    AsmMetric_measure(&c,&ref,rotated,&trial);
    if (!AsmMetric_preserved(&original,&trial)) fails++;
    uv[2] = -10; AsmMetric_measure(&c,&ref,uv,&trial);
    if (!trial.invalid || AsmMetric_preserved(&original,&trial)) fails++;
    /* A thin original triangle must not vanish in a Gram subtraction. */
    xyz[6] = 10; xyz[7] = 1e-6f; uv[2] = uv[4] = 10; uv[5] = 1e-6f;
    if (AsmMetric_prepare(&c,&ref)) fails++;
    AsmMetric_measure(&c,&ref,uv,&trial);
    if (trial.invalid || fabs(trial.minimum-1) > 1e-10 || fabs(trial.maximum-1) > 1e-10) fails++;
    {
        /* A tiny face can cross K while aggregate chart-area gates still pass.
         * Encoding must retain that original face's individual certificate. */
        float x[] = {0,0,0,10,0,0,0,10,0,20,0,0,21,0,0,20,1e-6f,0};
        float u[] = {0,0,10,0,0,10,20,0,21,0,20,1e-6f}; int32_t tri[] = {0,1,2,3,4,5};
        double encoded[12]; for (int k=0;k<12;k++) encoded[k]=u[k]; encoded[11]*=.74;
        AsmChart chart={0};chart.nv=6;chart.nf=2;chart.xyz=x;chart.uv=u;chart.faces=tri;AsmMetricFace metric[2];
        if (AsmMetric_prepare(&chart,metric)) fails++;
        AsmMetric_measure(&chart,metric,u,&original);AsmMetric_measure_encoded_uv64(&chart,metric,encoded,&trial);
        if (!AsmMetric_preserved(&original,&trial) || AsmMetric_encoded_faces_preserved(&chart,metric,encoded)) fails++;
        float local[12];for(int k=0;k<12;k++)local[k]=(float)encoded[k];
        if(!AsmMetric_local_face_preserved(&chart,0,metric,local) || AsmMetric_local_face_preserved(&chart,1,metric+1,local))fails++;
        /* Endpoint metrics and orientation agree for a half-turn, but its
         * linear deformation collapses at the midpoint. */
        for(int k=0;k<12;k++)local[k]=-u[k];
        if(AsmMetric_local_face_preserved(&chart,0,metric,local))fails++;
        for(int k=0;k<12;k++)local[k]=u[k];
        if(!AsmMetric_local_face_preserved(&chart,0,metric,local) || !AsmMetric_local_face_preserved(&chart,1,metric+1,local))fails++;
    }
    {
        /* A presentation reflection of placed coordinates is certified only
         * together with the parity it declares: u -> -u with the mirror bit
         * toggled keeps every face; the same coordinates under the old parity
         * are all flipped. */
        float x[] = {0,0,0, 10,0,0, 0,10,0}, u[] = {0,0, 10,0, 0,10}; int32_t tri[] = {0,1,2};
        double placed[6] = {100,5, 110,5, 100,15}, reflected[6];
        AsmChart chart={0};AsmMetricFace metric[1];
        chart.nv=3;chart.nf=1;chart.xyz=x;chart.uv=u;chart.faces=tri;chart.placed_uv=placed;
        for (int k=0;k<3;k++) { reflected[2*k] = 200-placed[2*k]; reflected[2*k+1] = placed[2*k+1]; }
        if (AsmMetric_prepare(&chart,metric)) fails++;
        if (!AsmMetric_encoded_faces_preserved(&chart,metric,placed)) fails++;
        if (AsmMetric_encoded_faces_preserved(&chart,metric,reflected)) fails++;
        if (!AsmMetric_encoded_faces_preserved_parity(&chart,metric,reflected,1)) fails++;
        if (AsmMetric_encoded_faces_preserved_parity(&chart,metric,placed,1)) fails++;
        chart.flags = ASM_CHART_MIRROR;          /* and back: a mirrored chart un-reflected to parity 0 */
        {
            double mirrored_placed[6]; for (int k=0;k<6;k++) mirrored_placed[k] = reflected[k];
            chart.placed_uv = mirrored_placed;
            if (!AsmMetric_encoded_faces_preserved(&chart,metric,mirrored_placed)) fails++;
            if (!AsmMetric_encoded_faces_preserved_parity(&chart,metric,placed,0)) fails++;
        }
    }
    fprintf(stderr,"  assembly original chart metric: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}
