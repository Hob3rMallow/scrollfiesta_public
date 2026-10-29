/* c3_stretch_metric.StretchBarrier / StretchPressure. Envelopes belong to
 * individual original moving faces, and tighten only between line searches. */
struct ArsStretch {double lower,upper,weight;int initialized,lower_qualified,upper_qualified;};
typedef struct ArsPressure {double weight,reference;int initialized,iterations;} ArsPressure;

static int ars_singular(const double j[4],double *lo,double *hi)
{
    for(int k=0;k<4;k++)if(!isfinite(j[k]))return 0;
    double scale=0;for(int k=0;k<4;k++)scale=fmax(scale,fabs(j[k]));if(!(scale>0))return 0;
    double a=j[0]/scale,b=j[1]/scale,c=j[2]/scale,d=j[3]/scale,det=a*d-b*c;
    if(!(det>0))return 0;
    double largest=.5*(hypot(a+d,c-b)+hypot(a-d,c+b));
    *hi=scale*largest;*lo=scale*(det/largest);return *lo>0 && isfinite(*hi);
}
static int ars_stretch_update(ArsStretch *band,const double j[4],double weight)
{
    double lo,hi;if(!ars_singular(j,&lo,&hi))return 0;
    double lower=lo>.75?.75:lo/1.05,upper=hi<1.25?1.25:hi*1.05;
    if(band->initialized){
        if((lo<.75 && band->lower_qualified) || (hi>1.25 && band->upper_qualified))return 0;
        lower=fmax(lower,band->lower);upper=fmin(upper,band->upper);
    }
    if(!(lo>lower && hi<upper) || !isfinite(upper))return 0;
    band->lower=lower;band->upper=upper;band->weight=weight;band->initialized=1;
    band->lower_qualified|=lo>=.75;band->upper_qualified|=hi<=1.25;return 1;
}
static double ars_stretch_face(const ArsStretch *band,const double j[4],double area,double gradient[4],double hessian[16])
{
    double lo,hi;if(!band->initialized || !ars_singular(j,&lo,&hi) ||
        !(lo>band->lower && hi<band->upper) ||
        (band->lower_qualified && lo<.75) || (band->upper_qualified && hi>1.25))return INFINITY;
    double a=1-band->lower,b=band->upper-1,c=8/(1/(a*a)+1/(b*b));
    double singular[]={hi,lo},first[2],second[2],energy=0,weight=area*band->weight;
    for(int k=0;k<2;k++){
        double t=singular[k]-1,dl=singular[k]-band->lower,du=band->upper-singular[k];
        energy+=c*(-log1p(t/a)-log1p(-t/b)+t*(1/a-1/b));
        first[k]=c*(-1/dl+1/du+1/a-1/b);second[k]=c*(1/(dl*dl)+1/(du*du));
    }
    if(gradient || hessian){
        double scale=0;for(int k=0;k<4;k++)scale=fmax(scale,fabs(j[k]));
        double x=j[0]/scale,y=j[1]/scale,z=j[2]/scale,w=j[3]/scale;
        double angle=.5*atan2(2*(x*y+z*w),x*x+z*z-y*y-w*w),vx=cos(angle),vy=sin(angle);
        double ux=x*vx+y*vy,uy=z*vx+w*vy,length=hypot(ux,uy);if(!(length>0))return INFINITY;
        ux/=length;uy/=length;
        double q[4][4]={{ux*vx,ux*vy,uy*vx,uy*vy},{uy*vy,-uy*vx,-ux*vy,ux*vx}};
        const double root2=1.4142135623730950488;
        double cross0[]={-ux*vy,ux*vx,-uy*vy,uy*vx},cross1[]={-uy*vx,-uy*vy,ux*vx,ux*vy};
        for(int k=0;k<4;k++){
            q[2][k]=(cross0[k]+cross1[k])/root2;q[3][k]=(cross0[k]-cross1[k])/root2;
            if(gradient)gradient[k]=weight*(first[0]*q[0][k]+first[1]*q[1][k]);
        }
        if(hessian){
            double eigen[]={second[0],second[1],c*(1/((hi-band->lower)*(lo-band->lower))+
                1/((band->upper-hi)*(band->upper-lo))),(first[0]+first[1])/(hi+lo)};
            memset(hessian,0,16*sizeof(double));
            for(int mode=0;mode<4;mode++)for(int r=0;r<4;r++)for(int col=0;col<4;col++)
                hessian[4*r+col]+=weight*fmax(eigen[mode],0)*q[mode][r]*q[mode][col];
        }
    }
    return weight*energy;
}
static int ars_pressure_increase(ArsPressure *p,double violation)
{
    if(!(violation>0) || p->weight>=1e8)return 0;
    p->weight=fmin(1e8,4*p->weight);p->reference=violation;return 1;
}
static int ars_pressure_begin(ArsPressure *p,double violation)
{
    int increased=0;
    if(!p->initialized){p->initialized=1;p->reference=violation;}
    else if(p->iterations%6==0){
        if(violation>0 && violation>.9*p->reference)increased=ars_pressure_increase(p,violation);
        p->reference=violation;
    }
    p->iterations++;return increased;
}
static int ars_stretch_begin(ArsProblem *p,ArsPressure *pressure,double *violation)
{
    AprSolver *s=p->s;AsmField *f=s->field;*violation=0;
    for(size_t at=0;at<s->n_active;at++){
        size_t t=s->active_faces[at];double j[4],lo,hi,det;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
        if(!ars_singular(j,&lo,&hi))return 0;*violation=fmax(*violation,fmax(.75/lo-1,hi/1.25-1));
    }
    if(ars_pressure_begin(pressure,*violation))s->stats->stretch_pressure_increases++;
    for(size_t at=0;at<s->n_active;at++){
        size_t t=s->active_faces[at];double j[4],lo,hi,det;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
        if(!ars_stretch_update(p->stretch+at,j,pressure->weight))return 0;
    }
    return 1;
}
