#include "shaft_warp.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double point[3], width, parameter, arc, length;
    /* C(q) = point + width*(b*q+c*q^2+d*q^3), q in [0,1]. */
    double coefficient[3][3];
    double lower[3], upper[3];
    double frame[5][3], quarter_arc[5];
} ShaftSpan;

struct ShaftWarp {
    size_t count;
    double *points;
    ShaftSpan *spans;
    double parameter_length, arc_length;
};

uint32_t ShaftWarp_abi(void) { return 1; }
size_t ShaftWarp_projection_size(void) { return sizeof(ShaftWarpProjection); }

static double sw_dot(const double a[3], const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static double sw_norm(const double a[3]) { return sqrt(sw_dot(a,a)); }
static void sw_cross(const double a[3], const double b[3], double out[3])
{
    out[0]=a[1]*b[2]-a[2]*b[1]; out[1]=a[2]*b[0]-a[0]*b[2];
    out[2]=a[0]*b[1]-a[1]*b[0];
}
static int sw_finite3(const double *a)
{ return a && isfinite(a[0]) && isfinite(a[1]) && isfinite(a[2]); }

static double sw_poly(const double *c, int degree, double q)
{
    double value=c[degree];
    while (degree>0) value=value*q+c[--degree];
    return value;
}
static int sw_root_add(double *roots, int *count, double root)
{
    if (root < -1e-12 || root > 1.+1e-12) return 0;
    root=fmin(1.,fmax(0.,root));
    for (int i=0;i<*count;i++) if (fabs(roots[i]-root)<2e-12) return 0;
    /* Near-multiple roots can also admit approximate critical points. Refuse
     * an unresolved set instead of assuming the exact polynomial root bound. */
    if (*count>=10) return -1;
    int at=*count;
    while (at>0 && roots[at-1]>root) { roots[at]=roots[at-1]; at--; }
    roots[at]=root; (*count)++; return 0;
}

/* Derivative roots partition the interval into monotone pieces. Testing the
 * critical points also retains even-multiplicity roots without sign changes. */
static int sw_roots(const double *input, int degree, double *roots)
{
    double c[10], derivative[10], critical[10], limits[12], scale=0.;
    int count=0, nc;
    for (int i=0;i<=degree;i++) scale=fmax(scale,fabs(input[i]));
    if (!(scale>0.) || !isfinite(scale)) return 0;
    for (int i=0;i<=degree;i++) c[i]=input[i]/scale;
    while (degree>0 && fabs(c[degree])<1e-14) degree--;
    if (degree==0) return 0;
    if (degree==1) { if (sw_root_add(roots,&count,-c[0]/c[1])) return -1; return count; }
    for (int i=1;i<=degree;i++) derivative[i-1]=i*c[i];
    nc=sw_roots(derivative,degree-1,critical);
    if (nc<0) return -1;
    limits[0]=0.;
    for (int i=0;i<nc;i++) limits[i+1]=critical[i];
    limits[nc+1]=1.;
    for (int i=0;i<nc+2;i++)
        if (fabs(sw_poly(c,degree,limits[i]))<=2e-13)
            if (sw_root_add(roots,&count,limits[i])) return -1;
    for (int i=0;i<nc+1;i++) {
        double lo=limits[i], hi=limits[i+1];
        double flo=sw_poly(c,degree,lo), fhi=sw_poly(c,degree,hi);
        if (!(hi>lo) || flo==0. || fhi==0. || (flo<0.)==(fhi<0.)) continue;
        for (int iteration=0;iteration<60 && hi-lo>2e-14;iteration++) {
            double mid=.5*(lo+hi), fm=sw_poly(c,degree,mid);
            if (fm==0.) { lo=hi=mid; break; }
            if ((fm<0.)==(flo<0.)) { lo=mid; flo=fm; }
            else hi=mid;
        }
        if (sw_root_add(roots,&count,.5*(lo+hi))) return -1;
    }
    return count;
}

static void sw_jet(const ShaftSpan *span, double q, double position[3],
                   double velocity[3], double acceleration_q[3])
{
    for (int a=0;a<3;a++) {
        double b=span->coefficient[a][0], c=span->coefficient[a][1], d=span->coefficient[a][2];
        if (position) position[a]=span->point[a]+span->width*q*(b+q*(c+q*d));
        velocity[a]=b+q*(2*c+3*q*d);
        if (acceleration_q) acceleration_q[a]=2*c+6*q*d;
    }
}
static double sw_speed(const ShaftSpan *span, double q)
{ double v[3]; sw_jet(span,q,NULL,v,NULL); return sw_norm(v); }

static double sw_gauss8(const ShaftSpan *span, double lo, double hi)
{
    static const double node[4]={.183434642495649805,.525532409916328986,
                                .796666477413626740,.960289856497536232};
    static const double weight[4]={.362683783378361983,.313706645877887287,
                                  .222381034453374471,.101228536290376259};
    double center=.5*(lo+hi), half=.5*(hi-lo), sum=0.;
    for (int i=0;i<4;i++) sum+=weight[i]*(sw_speed(span,center-half*node[i])+sw_speed(span,center+half*node[i]));
    return half*sum;
}
static int sw_integral(const ShaftSpan *span, double lo, double hi, double tolerance,
                       int depth, double *out)
{
    if (lo==hi) { *out=0.; return 0; }
    double mid=.5*(lo+hi), whole=sw_gauss8(span,lo,hi);
    double left=sw_gauss8(span,lo,mid), right=sw_gauss8(span,mid,hi);
    if (!isfinite(whole+left+right)) return -3;
    if (fabs(whole-left-right)<=fmax(tolerance,2e-15*(hi-lo))) { *out=left+right; return 0; }
    if (depth==0 || mid==lo || mid==hi) return -3;
    if (sw_integral(span,lo,mid,tolerance/2,depth-1,&left) ||
        sw_integral(span,mid,hi,tolerance/2,depth-1,&right)) return -3;
    *out=left+right; return 0;
}

static int sw_frame_rhs(const ShaftSpan *span, double q, const double normal[3], double out[3])
{
    double v[3], a[3], tangent[3], derivative[3], speed, along;
    sw_jet(span,q,NULL,v,a); speed=sw_norm(v);
    if (!(speed>0.) || !isfinite(speed)) return -3;
    for (int i=0;i<3;i++) tangent[i]=v[i]/speed;
    along=sw_dot(tangent,a);
    for (int i=0;i<3;i++) derivative[i]=(a[i]-tangent[i]*along)/speed;
    along=sw_dot(normal,derivative);
    for (int i=0;i<3;i++) out[i]=-tangent[i]*along;
    return 0;
}
static int sw_frame_step(const ShaftSpan *span, double lo, double hi, const double start[3], double out[3])
{
    double k1[3],k2[3],k3[3],k4[3],work[3],width=hi-lo;
    if (sw_frame_rhs(span,lo,start,k1)) return -3;
    for (int a=0;a<3;a++) work[a]=start[a]+width*.5*k1[a];
    if (sw_frame_rhs(span,.5*(lo+hi),work,k2)) return -3;
    for (int a=0;a<3;a++) work[a]=start[a]+width*.5*k2[a];
    if (sw_frame_rhs(span,.5*(lo+hi),work,k3)) return -3;
    for (int a=0;a<3;a++) work[a]=start[a]+width*k3[a];
    if (sw_frame_rhs(span,hi,work,k4)) return -3;
    for (int a=0;a<3;a++) out[a]=start[a]+width*(k1[a]+2*k2[a]+2*k3[a]+k4[a])/6;
    return 0;
}
static int sw_frame_clean(const ShaftSpan *span, double q, double normal[3])
{
    double v[3],speed,projection,norm;
    sw_jet(span,q,NULL,v,NULL); speed=sw_norm(v);
    if (!(speed>0.) || !isfinite(speed)) return -3;
    projection=sw_dot(v,normal)/(speed*speed);
    for (int a=0;a<3;a++) normal[a]-=projection*v[a];
    norm=sw_norm(normal);
    if (!(norm>1e-12) || !isfinite(norm)) return -3;
    for (int a=0;a<3;a++) normal[a]/=norm;
    return 0;
}
static int sw_frame_integral(const ShaftSpan *span, double lo, double hi, const double start[3],
                             double tolerance, int depth, double out[3])
{
    double whole[3],half[3],fine[3],difference[3],mid=.5*(lo+hi);
    if (lo==hi) { memcpy(out,start,3*sizeof(double)); return 0; }
    if (sw_frame_step(span,lo,hi,start,whole) || sw_frame_step(span,lo,mid,start,half) ||
        sw_frame_step(span,mid,hi,half,fine)) return -3;
    for (int a=0;a<3;a++) difference[a]=whole[a]-fine[a];
    if (sw_norm(difference)<=fmax(tolerance,2e-14)) {
        memcpy(out,fine,3*sizeof(double)); return sw_frame_clean(span,hi,out);
    }
    if (depth==0 || mid==lo || mid==hi) return -3;
    if (sw_frame_integral(span,lo,mid,start,tolerance/2,depth-1,half)) return -3;
    return sw_frame_integral(span,mid,hi,half,tolerance/2,depth-1,out);
}

void ShaftWarp_free(ShaftWarp *warp)
{
    if (warp) { free(warp->points); free(warp->spans); free(warp); }
}

int ShaftWarp_create(ShaftWarp **out, const double *points, size_t count, const double initial_normal[3])
{
    ShaftWarp *warp=NULL;
    double *width=NULL,*diagonal=NULL,*rhs=NULL,*second=NULL;
    int status=-2;
    if (!out || *out || !points || count<3 || count>SIZE_MAX/(3*sizeof(double)) ||
        count>SIZE_MAX/sizeof(ShaftSpan) || (initial_normal && !sw_finite3(initial_normal))) return -1;
    for (size_t i=0;i<count;i++) if (!sw_finite3(points+3*i)) return -1;
    warp=(ShaftWarp*)calloc(1,sizeof(*warp));
    if (!warp) return -2;
    warp->count=count;
    warp->points=(double*)malloc(3*count*sizeof(double));
    warp->spans=(ShaftSpan*)calloc(count-1,sizeof(ShaftSpan));
    width=(double*)malloc((count-1)*sizeof(double));
    diagonal=(double*)calloc(count,sizeof(double));
    rhs=(double*)calloc(3*count,sizeof(double));
    second=(double*)calloc(3*count,sizeof(double));
    if (!warp->points || !warp->spans || !width || !diagonal || !rhs || !second) goto done;
    memcpy(warp->points,points,3*count*sizeof(double));
    status=-1;
    for (size_t i=0;i<count-1;i++) {
        double delta[3]; for (int a=0;a<3;a++) delta[a]=points[3*(i+1)+a]-points[3*i+a];
        width[i]=sw_norm(delta);
        if (!(width[i]>1e-10) || !isfinite(width[i])) goto done;
    }
    for (size_t i=1;i<count-1;i++) {
        diagonal[i]=2*(width[i-1]+width[i]);
        for (int a=0;a<3;a++) rhs[3*i+a]=6*((points[3*(i+1)+a]-points[3*i+a])/width[i]
                                                  -(points[3*i+a]-points[3*(i-1)+a])/width[i-1]);
        if (i>1) {
            double factor=width[i-1]/diagonal[i-1];
            diagonal[i]-=factor*width[i-1];
            for (int a=0;a<3;a++) rhs[3*i+a]-=factor*rhs[3*(i-1)+a];
        }
    }
    for (size_t i=count-2;i>0;i--)
        for (int a=0;a<3;a++) second[3*i+a]=(rhs[3*i+a]-width[i]*second[3*(i+1)+a])/diagonal[i];
    status=-3;
    for (size_t i=0;i<count-1;i++) {
        ShaftSpan *span=warp->spans+i;
        double speed2[5]={0}, derivative[4], roots[10], checks[12], first[3];
        int nr;
        span->width=width[i]; span->parameter=warp->parameter_length; span->arc=warp->arc_length;
        for (int a=0;a<3;a++) {
            double m=second[3*i+a], next=second[3*(i+1)+a], *c=span->coefficient[a];
            span->point[a]=points[3*i+a];
            c[0]=(points[3*(i+1)+a]-points[3*i+a])/width[i]-width[i]*(2*m+next)/6;
            c[1]=width[i]*m/2; c[2]=width[i]*(next-m)/6;
            double v[3]={c[0],2*c[1],3*c[2]};
            for (int j=0;j<3;j++) for (int k=0;k<3;k++) speed2[j+k]+=v[j]*v[k];
            double control[4]={points[3*i+a],points[3*i+a]+width[i]*c[0]/3,
                points[3*(i+1)+a]-width[i]*(c[0]+2*c[1]+3*c[2])/3,points[3*(i+1)+a]};
            span->lower[a]=span->upper[a]=control[0];
            for (int j=1;j<4;j++) { span->lower[a]=fmin(span->lower[a],control[j]); span->upper[a]=fmax(span->upper[a],control[j]); }
        }
        for (int j=1;j<5;j++) derivative[j-1]=j*speed2[j];
        nr=sw_roots(derivative,3,roots);
        if (nr<0) goto done;
        checks[0]=0.; checks[nr+1]=1.;
        for (int j=0;j<nr;j++) checks[j+1]=roots[j];
        for (int j=0;j<nr+2;j++) if (!(sw_speed(span,checks[j])>1e-8)) goto done;
        if (i>0) memcpy(first,warp->spans[i-1].frame[4],sizeof(first));
        else if (initial_normal) memcpy(first,initial_normal,sizeof(first));
        else {
            double v[3]; sw_jet(span,0.,NULL,v,NULL);
            first[0]=0.; first[1]=1.; first[2]=0.;
            if (fabs(v[1])/sw_norm(v)>.95) {
                int least=0; for (int a=1;a<3;a++) if (fabs(v[a])<fabs(v[least])) least=a;
                first[0]=first[1]=first[2]=0.; first[least]=1.;
            }
        }
        if (sw_frame_clean(span,0.,first)) goto done;
        memcpy(span->frame[0],first,sizeof(first));
        for (int quarter=1;quarter<=4;quarter++) {
            double integral, lo=(quarter-1)*.25, hi=quarter*.25;
            if (sw_integral(span,lo,hi,2.5e-13,32,&integral) ||
                sw_frame_integral(span,lo,hi,span->frame[quarter-1],2.5e-12,32,span->frame[quarter])) goto done;
            span->quarter_arc[quarter]=span->quarter_arc[quarter-1]+width[i]*integral;
        }
        span->length=span->quarter_arc[4];
        warp->parameter_length+=width[i]; warp->arc_length+=span->length;
        if (!isfinite(warp->arc_length+warp->parameter_length)) goto done;
    }
    *out=warp; warp=NULL; status=0;
done:
    free(width); free(diagonal); free(rhs); free(second); ShaftWarp_free(warp);
    return status;
}

double ShaftWarp_length_um(const ShaftWarp *warp) { return warp ? warp->arc_length : NAN; }
double ShaftWarp_parameter_length_um(const ShaftWarp *warp) { return warp ? warp->parameter_length : NAN; }

static int sw_arc_local(const ShaftSpan *span, double q, double *out)
{
    int quarter=(int)floor(q*4); double integral;
    if (quarter>=4) { *out=span->length; return 0; }
    if (sw_integral(span,quarter*.25,q,2.5e-13,32,&integral)) return -3;
    *out=span->quarter_arc[quarter]+span->width*integral; return 0;
}
static int sw_evaluate(const ShaftWarp *warp, size_t index, double q, ShaftWarpProjection *out)
{
    const ShaftSpan *span=warp->spans+index;
    ShaftWarpProjection value;
    double v[3],a[3],speed,along,arc;
    int quarter=(int)floor(q*4);
    memset(&value,0,sizeof(value));
    sw_jet(span,q,value.foot_um_zyx,v,a); speed=sw_norm(v);
    if (!(speed>0.) || !isfinite(speed)) return -3;
    for (int k=0;k<3;k++) value.tangent_zyx[k]=v[k]/speed;
    along=sw_dot(value.tangent_zyx,a);
    for (int k=0;k<3;k++) value.curvature_per_um_zyx[k]=(a[k]-along*value.tangent_zyx[k])/(span->width*speed*speed);
    if (quarter>=4) memcpy(value.normal1_zyx,span->frame[4],3*sizeof(double));
    else if (sw_frame_integral(span,quarter*.25,q,span->frame[quarter],2.5e-12,32,value.normal1_zyx)) return -3;
    sw_cross(value.tangent_zyx,value.normal1_zyx,value.normal2_zyx);
    if (sw_arc_local(span,q,&arc)) return -3;
    value.parameter_um=span->parameter+q*span->width; value.s_um=span->arc+arc;
    value.local_jacobian=1.;
    if (value.parameter_um<=1e-8 || value.parameter_um>=warp->parameter_length-1e-8) value.flags|=SHAFT_WARP_ENDPOINT;
    *out=value; return 0;
}

int ShaftWarp_eval_parameter(const ShaftWarp *warp, double parameter, ShaftWarpProjection *out)
{
    size_t lo=0,hi;
    if (!warp || !out || !isfinite(parameter) || parameter<0. || parameter>warp->parameter_length) return -1;
    hi=warp->count-1;
    while (hi-lo>1) { size_t mid=lo+(hi-lo)/2; if (warp->spans[mid].parameter<=parameter) lo=mid; else hi=mid; }
    const ShaftSpan *span=warp->spans+lo;
    return sw_evaluate(warp,lo,fmin(1.,fmax(0.,(parameter-span->parameter)/span->width)),out);
}

int ShaftWarp_eval_arc(const ShaftWarp *warp, double arc, ShaftWarpProjection *out)
{
    size_t lo=0,hi;
    double a=0.,b=1.,q,target;
    if (!warp || !out || !isfinite(arc) || arc<0. || arc>warp->arc_length) return -1;
    hi=warp->count-1;
    while (hi-lo>1) { size_t mid=lo+(hi-lo)/2; if (warp->spans[mid].arc<=arc) lo=mid; else hi=mid; }
    const ShaftSpan *span=warp->spans+lo;
    target=arc-span->arc;
    if (target<=1e-9) return sw_evaluate(warp,lo,0.,out);
    if (span->length-target<=1e-9) return sw_evaluate(warp,lo,1.,out);
    q=target/span->length;
    for (int iteration=0;iteration<64;iteration++) {
        double value,next;
        if (sw_arc_local(span,q,&value)) return -3;
        double error=value-target;
        if (fabs(error)<=1e-8 || b-a<=2e-14) return sw_evaluate(warp,lo,q,out);
        if (error<0.) a=q; else b=q;
        next=q-error/(span->width*sw_speed(span,q));
        q=(next>a && next<b) ? next : .5*(a+b);
    }
    return -3;
}

int ShaftWarp_bounds_arc(const ShaftWarp *warp, double begin, double end,
    double radius, double lower[3], double upper[3])
{
    ShaftWarpProjection first,last;
    double lo[3]={DBL_MAX,DBL_MAX,DBL_MAX},hi[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
    if (!warp || !lower || !upper || lower==upper ||
        !isfinite(begin) || !isfinite(end) || !isfinite(radius) ||
        begin<0. || end<begin || end>warp->arc_length || radius<0.) return -1;
    if (ShaftWarp_eval_arc(warp,begin,&first) || ShaftWarp_eval_arc(warp,end,&last)) return -3;
    for (size_t i=0;i<warp->count-1;i++) {
        const ShaftSpan *span=warp->spans+i;
        if (span->parameter>last.parameter_um ||
            span->parameter+span->width<first.parameter_um) continue;
        double a=fmax(0.,(first.parameter_um-span->parameter)/span->width);
        double b=fmin(1.,(last.parameter_um-span->parameter)/span->width);
        if (a>b) continue;
        for (int axis=0;axis<3;axis++) {
            const double *c=span->coefficient[axis];
            double velocity[3]={c[0],2*c[1],3*c[2]},roots[10],queries[12]={a,b};
            int nr=sw_roots(velocity,2,roots),nq=2;
            if (nr<0) return -3;
            for (int k=0;k<nr;k++) if (roots[k]>a && roots[k]<b) queries[nq++]=roots[k];
            for (int k=0;k<nq;k++) {
                double q=queries[k],value=span->point[axis]+span->width*q*(c[0]+q*(c[1]+q*c[2]));
                if (!isfinite(value)) return -3;
                lo[axis]=fmin(lo[axis],value); hi[axis]=fmax(hi[axis],value);
            }
        }
    }
    for (int axis=0;axis<3;axis++) {
        /* Cover arc-inversion tolerance and coordinate arithmetic roundoff.
         * A disk's displacement along any world axis is at most its radius. */
        double slack=1e-7+128*DBL_EPSILON*(fmax(fabs(lo[axis]),fabs(hi[axis]))+
                                           warp->parameter_length+radius);
        if (lo[axis]>hi[axis]) return -3;
        lo[axis]-=radius+slack; hi[axis]+=radius+slack;
        if (!isfinite(lo[axis]) || !isfinite(hi[axis])) return -3;
    }
    memcpy(lower,lo,sizeof lo); memcpy(upper,hi,sizeof hi); return 0;
}

typedef struct { size_t span; double lower; } ShaftSearch;
typedef struct { size_t span; double q,distance; } ShaftFoot;
static int sw_search_order(const void *a, const void *b)
{
    const ShaftSearch *left=(const ShaftSearch*)a,*right=(const ShaftSearch*)b;
    return left->lower<right->lower ? -1 : left->lower>right->lower ? 1 : left->span<right->span ? -1 : left->span>right->span;
}

int ShaftWarp_project(const ShaftWarp *warp, const double point[3], double ambiguity, ShaftWarpProjection *out)
{
    ShaftSearch *search=NULL;
    ShaftFoot *feet=NULL;
    size_t nfeet=0,best=0,nspans;
    double upper=DBL_MAX,best_distance=DBL_MAX;
    int status=-2;
    if (!warp || !out || !sw_finite3(point) || !isfinite(ambiguity) || ambiguity<0.) return -1;
    nspans=warp->count-1;
    if (nspans>(SIZE_MAX/sizeof(ShaftFoot)-2)/10) return -2;
    search=(ShaftSearch*)malloc(nspans*sizeof(*search));
    feet=(ShaftFoot*)malloc((10*nspans+2)*sizeof(*feet));
    if (!search || !feet) goto done;
    for (size_t i=0;i<warp->count;i++) {
        double delta[3]; for (int a=0;a<3;a++) delta[a]=warp->points[3*i+a]-point[a];
        upper=fmin(upper,sw_norm(delta));
    }
    for (size_t i=0;i<nspans;i++) {
        const ShaftSpan *span=warp->spans+i; double delta[3];
        for (int a=0;a<3;a++) delta[a]=fmax(0.,fmax(span->lower[a]-point[a],point[a]-span->upper[a]));
        search[i].span=i; search[i].lower=sw_norm(delta);
    }
    qsort(search,nspans,sizeof(*search),sw_search_order);
    for (size_t order=0;order<nspans;order++) {
        size_t i=search[order].span; const ShaftSpan *span=warp->spans+i;
        double polynomial[6]={0},derivative[5],roots[10],queries[12],magnitude=0.;
        int nr,nq=0;
        if (search[order].lower>upper+ambiguity+1e-7) break;
        for (int a=0;a<3;a++) {
            const double *c=span->coefficient[a];
            double position[4]={(span->point[a]-point[a])/span->width,c[0],c[1],c[2]};
            double velocity[3]={c[0],2*c[1],3*c[2]};
            for (int j=0;j<4;j++) for (int k=0;k<3;k++) polynomial[j+k]+=position[j]*velocity[k];
        }
        for (int k=0;k<6;k++) magnitude=fmax(magnitude,fabs(polynomial[k]));
        for (int k=1;k<6;k++) derivative[k-1]=k*polynomial[k];
        nr=sw_roots(polynomial,5,roots);
        if (nr<0) { status=-3; goto done; }
        if (i==0 && polynomial[0]>=0.) queries[nq++]=0.;
        for (int k=0;k<nr;k++) if (sw_poly(derivative,4,roots[k])>=-1e-10*fmax(1.,magnitude)) queries[nq++]=roots[k];
        if (i==nspans-1 && sw_poly(polynomial,5,1.)<=0.) queries[nq++]=1.;
        for (int k=0;k<nq;k++) {
            double position[3],velocity[3],delta[3];
            sw_jet(span,queries[k],position,velocity,NULL);
            for (int a=0;a<3;a++) delta[a]=position[a]-point[a];
            double distance=sw_norm(delta);
            if (!isfinite(distance)) { status=-3; goto done; }
            feet[nfeet].span=i; feet[nfeet].q=queries[k]; feet[nfeet].distance=distance;
            if (distance<best_distance) { best_distance=distance; best=nfeet; }
            nfeet++; upper=fmin(upper,distance);
        }
    }
    status=-3;
    if (!nfeet) goto done;
    ShaftWarpProjection value;
    if (sw_evaluate(warp,feet[best].span,feet[best].q,&value)) goto done;
    double min_arc=value.s_um,max_arc=value.s_um,offset[3];
    for (size_t i=0;i<nfeet;i++) if (feet[i].distance<=best_distance+ambiguity) {
        const ShaftSpan *span=warp->spans+feet[i].span; double arc;
        if (sw_arc_local(span,feet[i].q,&arc)) goto done;
        arc+=span->arc; min_arc=fmin(min_arc,arc); max_arc=fmax(max_arc,arc);
    }
    if (max_arc-min_arc>fmax(1e-6,2*ambiguity)) value.flags|=SHAFT_WARP_AMBIGUOUS;
    for (int a=0;a<3;a++) offset[a]=point[a]-value.foot_um_zyx[a];
    value.distance_um=best_distance; value.normal_residual_um=sw_dot(offset,value.tangent_zyx);
    value.local_jacobian=1.-sw_dot(offset,value.curvature_per_um_zyx);
    if (value.local_jacobian<=0.) value.flags|=SHAFT_WARP_FOLDED;
    if ((value.flags&SHAFT_WARP_ENDPOINT) && fabs(value.normal_residual_um)>fmax(1e-6,64*DBL_EPSILON*warp->arc_length))
        value.flags|=SHAFT_WARP_OUTSIDE_SUPPORT;
    *out=value; status=0;
done:
    free(search); free(feet); return status;
}

static int sw_controls(double radius,double jacobian,double ambiguity)
{ return isfinite(radius) && radius>=0. && isfinite(jacobian) && jacobian>0. && jacobian<=1. && isfinite(ambiguity) && ambiguity>=0.; }
static uint32_t sw_domain(const ShaftWarpProjection *p,double radius,double jacobian)
{
    uint32_t flags=p->flags;
    /* Polynomial projection may move a radius-boundary foot by a few ulps.
     * This is numerical slack in um, not an uncertainty or sampling margin. */
    if (p->distance_um>radius+fmax(1e-7,64*DBL_EPSILON*radius)) flags|=SHAFT_WARP_OUTSIDE_RADIUS;
    if (p->local_jacobian<jacobian) flags|=SHAFT_WARP_BELOW_JACOBIAN;
    return flags;
}

int ShaftWarp_to_metric(const ShaftWarp *warp,const double point[3],double radius,double jacobian,double ambiguity,
                        double out[3],uint32_t *flags)
{
    ShaftWarpProjection p; double offset[3],value[3]; int status;
    if (!out || !flags || !sw_controls(radius,jacobian,ambiguity)) return -1;
    status=ShaftWarp_project(warp,point,ambiguity,&p); if (status) return status;
    *flags=sw_domain(&p,radius,jacobian);
    if (*flags & ~(uint32_t)SHAFT_WARP_ENDPOINT) return 1;
    for (int a=0;a<3;a++) offset[a]=point[a]-p.foot_um_zyx[a];
    value[0]=p.s_um; value[1]=sw_dot(offset,p.normal1_zyx); value[2]=sw_dot(offset,p.normal2_zyx);
    memcpy(out,value,sizeof(value)); return 0;
}

int ShaftWarp_from_metric(const ShaftWarp *warp,const double suv[3],double radius,double jacobian,double ambiguity,
                          double out[3],uint32_t *flags)
{
    ShaftWarpProjection frame,projected; double point[3]; int status;
    if (!warp || !out || !flags || !sw_finite3(suv) || !sw_controls(radius,jacobian,ambiguity)) return -1;
    *flags=0;
    if (suv[0]<0. || suv[0]>warp->arc_length) { *flags=SHAFT_WARP_OUTSIDE_SUPPORT; return 1; }
    status=ShaftWarp_eval_arc(warp,suv[0],&frame); if (status) return status;
    frame.distance_um=hypot(suv[1],suv[2]);
    for (int a=0;a<3;a++) {
        double offset=suv[1]*frame.normal1_zyx[a]+suv[2]*frame.normal2_zyx[a];
        point[a]=frame.foot_um_zyx[a]+offset;
        frame.local_jacobian-=offset*frame.curvature_per_um_zyx[a];
    }
    if (frame.local_jacobian<=0.) frame.flags|=SHAFT_WARP_FOLDED;
    *flags=sw_domain(&frame,radius,jacobian);
    if (*flags & ~(uint32_t)SHAFT_WARP_ENDPOINT) return 1;
    if (!sw_finite3(point) || !isfinite(frame.local_jacobian)) return -3;
    status=ShaftWarp_project(warp,point,ambiguity,&projected); if (status) return status;
    *flags=sw_domain(&projected,radius,jacobian);
    if (fabs(projected.s_um-suv[0])>fmax(1e-5,2*ambiguity)) *flags|=SHAFT_WARP_DIFFERENT_BRANCH;
    if (*flags & ~(uint32_t)SHAFT_WARP_ENDPOINT) return 1;
    memcpy(out,point,sizeof(point)); return 0;
}
