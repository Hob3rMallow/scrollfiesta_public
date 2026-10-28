#include "metric_bounds.h"
/* Private local metric constraints for the shared-source quadratic. These
 * tangent planes are not a nonlinear feasibility certificate: the unchanged
 * exact objective, face/path, physical and material gates still own acceptance.
 * Only originally good faces within 1e-4 of a 25% scale bound participate.
 * Their linear endpoint is directed 1e-4 inside that bound to avoid a tangent
 * direction whose second-order curvature would immediately cross the limit. */
typedef struct SmRest { double gx[3],gy[3],area; } SmRest;
typedef struct SmBound {
    int index[6]; double g[6],rhs; size_t face; int upper;
} SmBound;
typedef struct SmCurvature SmCurvature;
typedef struct SmResponseReuse SmResponseReuse;
typedef struct SmActiveReuse SmActiveReuse;
typedef struct SmLocalWarm SmLocalWarm;
typedef struct SmMetric {
    SmLocalWarm *ipm;
    const float *exact_xyz;
    size_t nf,constrained_calls,rows,active_peak,adds,drops,solves;
    const int32_t *faces; SmRest *rest; uint8_t *good,*feedback;
    const char *capture_prefix;
    int origin_policy;
    int equilibrate;double schur_diag_min,schur_diag_max;
    double max_violation;SmCurvature *curvature;SmResponseReuse *responses;SmActiveReuse *active_reuse;
    double profile_response,profile_schur_build,profile_schur_solve,profile_scan,profile_reconstruct;
    size_t profile_schur_calls,profile_gather_coordinates,profile_full_coordinates;
} SmMetric;

static int smm_exact(const SmMetric *m,size_t f,const double *uv,unsigned *bits)
{
    if(!m->exact_xyz||f>=m->nf)return -1;
    float xyz[9];double target[6];MetricBoundsResult result;
    for(int corner=0;corner<3;corner++) {
        size_t vertex=(size_t)m->faces[3*f+corner];
        memcpy(xyz+3*corner,m->exact_xyz+3*vertex,3*sizeof(float));
        memcpy(target+2*corner,uv+2*vertex,2*sizeof(double));
    }
    const double bounds[4]={1/1.1,1.1,1/1.25,1.25};
    if(MetricBounds_classify(xyz,target,bounds,&result))return -1;
    *bits=result.violated;return 0;
}

static int smm_singular(const double j[4],double sigma[2],double gradient[8])
{
    for(int k=0;k<4;k++)if(!isfinite(j[k]))return -1;
    double a=j[0]*j[0]+j[1]*j[1],b=j[0]*j[2]+j[1]*j[3],c=j[2]*j[2]+j[3]*j[3];
    double det=j[0]*j[3]-j[1]*j[2],theta=.5*atan2(2*b,a-c);
    sigma[1]=sqrt(fmax(0,.5*(a+c+hypot(a-c,2*b))));
    sigma[0]=sigma[1]>0?det/sigma[1]:0;
    if(!(sigma[0]>0)||!isfinite(sigma[1]))return -1;
    if(gradient)for(int k=0;k<2;k++) {
        double u=k?cos(theta):-sin(theta),v=k?sin(theta):cos(theta);
        double x=(j[0]*u+j[2]*v)/sigma[k],y=(j[1]*u+j[3]*v)/sigma[k];
        gradient[4*k]=u*x;gradient[4*k+1]=u*y;gradient[4*k+2]=v*x;gradient[4*k+3]=v*y;
    }
    return 0;
}
static void smm_jacobian(const SmRest *r,const int32_t *face,const double *uv,double j[4])
{
    memset(j,0,4*sizeof *j);
    for(size_t k=0;k<3;k++) {
        size_t v=(size_t)face[k];
        j[0]+=r->gx[k]*uv[2*v];j[1]+=r->gy[k]*uv[2*v];
        j[2]+=r->gx[k]*uv[2*v+1];j[3]+=r->gy[k]*uv[2*v+1];
    }
}
static int smm_setup(Arena_T arena,const float *xyz,size_t nv,const int32_t *faces,
                     size_t nf,const double *reference,SmMetric *m)
{
    if(!arena||!xyz||!nv||!faces||!nf||!reference||!m||nf>SIZE_MAX/2)return -1;
    memset(m,0,sizeof *m);m->nf=nf;m->faces=faces;
    m->rest=ARENA_ALLOC(arena,nf*sizeof *m->rest);m->good=ARENA_CALLOC(arena,nf,1);
    for(size_t f=0;f<nf;f++) {
        const int32_t *face=faces+3*f;double e[3],q[3],cross[3],length=0,projection=0;
        for(int k=0;k<3;k++)if(face[k]<0||(size_t)face[k]>=nv)return -1;
        for(int k=0;k<3;k++) {
            e[k]=(double)xyz[3*(size_t)face[1]+k]-xyz[3*(size_t)face[0]+k];
            q[k]=(double)xyz[3*(size_t)face[2]+k]-xyz[3*(size_t)face[0]+k];
            length+=e[k]*e[k];
        }
        length=sqrt(length);if(!(length>1e-12)||!isfinite(length))return -1;
        for(int k=0;k<3;k++)projection+=q[k]*(e[k]/length);
        cross[0]=e[1]*q[2]-e[2]*q[1];cross[1]=e[2]*q[0]-e[0]*q[2];cross[2]=e[0]*q[1]-e[1]*q[0];
        double height=hypot(hypot(cross[0],cross[1]),cross[2])/length;
        if(!(height>1e-12)||!isfinite(height))return -1;
        SmRest *r=m->rest+f;r->area=.5*length*height;r->gx[0]=-1/length;r->gx[1]=1/length;r->gx[2]=0;
        r->gy[0]=projection/(length*height)-1/height;r->gy[1]=-projection/(length*height);r->gy[2]=1/height;
        double j[4],sigma[2];smm_jacobian(r,face,reference,j);if(smm_singular(j,sigma,NULL))return -1;
        m->good[f]=(uint8_t)(sigma[0]>=.8&&sigma[1]<=1.25);
    }
    return 0;
}
static int smm_bounds(Arena_T arena,const SmMetric *m,const int32_t *map,
                      const double *uv,SmBound **out,size_t *count)
{
    SmBound *rows=ARENA_ALLOC(arena,2*m->nf*sizeof *rows);size_t n=0;
    for(size_t f=0;f<m->nf;f++)if(m->good[f]) {
        double j[4],sigma[2],gradient[8];const int32_t *face=m->faces+3*f;
        smm_jacobian(m->rest+f,face,uv,j);if(smm_singular(j,sigma,gradient))return -1;
        for(int side=0;side<2;side++)if((m->feedback&&(m->feedback[f]&(1u<<side)))||(side?sigma[1]>=1.25-1e-4:sigma[0]<=.8+1e-4)) {
            /* An arbitrary singular vector at a repeated singular value is
             * not the derivative of the minimum/maximum branch. Refuse it. */
            if(sigma[1]-sigma[0]<1e-10*fmax(1,sigma[1]))return -1;
            SmBound *r=rows+n++;r->face=f;r->upper=side;double sign=side?-1:1;
            r->rhs=side?sigma[1]-(1.25-1e-4):(.8+1e-4)-sigma[0];
            for(size_t k=0;k<3;k++)for(size_t xy=0;xy<2;xy++) {
                r->index[2*k+xy]=2*map[face[k]]+(int)xy;
                r->g[2*k+xy]=sign*(gradient[4*side+2*xy]*m->rest[f].gx[k]+
                                                   gradient[4*side+2*xy+1]*m->rest[f].gy[k]);
            }
        }
    }
    *out=rows;*count=n;return 0;
}
static double smm_dot(const SmBound *r,const double *x)
{
    double value=0;for(size_t k=0;k<6;k++)value+=r->g[k]*x[r->index[k]];return value;
}

#include "closeup_lg_origin.inc"
