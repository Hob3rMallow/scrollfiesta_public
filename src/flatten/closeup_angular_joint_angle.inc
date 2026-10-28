/* PRIVATE DIAGNOSTIC ONLY. The production objective is unchanged.
 * This touching Gauss--Newton model is not a global majorizer. Every proposal
 * must still pass the exact nonlinear objective and original geometric gates.
 * A joint complementary-axis assignment retains both raw physical arms. */
static double sra_axis_angle(const double p[2],int axis)
{
    double x=p[0],y=p[1];
    return .5*atan2(axis?-2*x*y:2*x*y,axis?y*y-x*x:x*x-y*y);
}

static int sra_residuals(const double h[2],const double v[2],double out[2])
{
    double a[2]={sra_axis_angle(h,0),sra_axis_angle(v,1)};
    double b[2]={sra_axis_angle(h,1),sra_axis_angle(v,0)};
    double aa=a[0]*a[0]+a[1]*a[1],bb=b[0]*b[0]+b[1]*b[1];
    double tolerance=16*DBL_EPSILON*fmax(1,fmax(aa,bb));
    int swapped=aa<bb-tolerance?0:bb<aa-tolerance?1:
        atan2(2*h[0]*h[1],h[0]*h[0]-h[1]*h[1])>
        atan2(2*v[0]*v[1],v[0]*v[0]-v[1]*v[1]);
    const double *selected=swapped?b:a;
    out[0]=selected[0];out[1]=selected[1];return swapped;
}

static double sra_energy(const SRDirection *d,const int32_t *face,const double *uv)
{
    if(!(d->weight>0))return 0;
    double h[2],v[2];sr_direction_vectors(d,face,uv,h,v);
    double hh=h[0]*h[0]+h[1]*h[1],vv=v[0]*v[0]+v[1]*v[1];
    if(!(hh>DBL_MIN && vv>DBL_MIN) || !isfinite(hh) || !isfinite(vv))return INFINITY;
    double a=sra_axis_angle(h,0),b=sra_axis_angle(v,1);
    double c=sra_axis_angle(h,1),e=sra_axis_angle(v,0);
    return d->weight*fmin(a*a+b*b,c*c+e*e);
}

/* Gradients are in local vertex coordinate order [u0,v0,u1,v1,u2,v2].
 * d theta / d p = (-p.y,p.x)/|p|^2 for either complementary output axis.
 * No direction is assigned to a missing-evidence face. */
static int sra_model(const SRDirection *d,const int32_t *face,const double *uv,
    double residual[2],double gradient[2][6])
{
    memset(residual,0,2*sizeof *residual);memset(gradient,0,12*sizeof(double));
    if(!(d->weight>0))return 0;
    double h[2],v[2];sr_direction_vectors(d,face,uv,h,v);
    double hh=h[0]*h[0]+h[1]*h[1],vv=v[0]*v[0]+v[1]*v[1];
    if(!(hh>DBL_MIN && vv>DBL_MIN) || !isfinite(hh) || !isfinite(vv))return -1;
    sra_residuals(h,v,residual);
    double normal[2][2]={{-h[1]/hh,h[0]/hh},{-v[1]/vv,v[0]/vv}};
    for(int arm=0;arm<2;arm++) {
        const double *coefficient=arm?d->v:d->h;
        for(int k=0;k<3;k++)for(int axis=0;axis<2;axis++)
            gradient[arm][2*k+axis]=coefficient[k]*normal[arm][axis];
    }
    return 0;
}

static int sra_assemble_face(SRSystem *system,size_t f,const double *uv)
{
    const SRDirection *d=system->direction+f;const int32_t *face=system->faces+3*f;
    double residual[2],g[2][6];if(sra_model(d,face,uv,residual,g))return -1;
    double intercept[2]={0,0};
    for(int arm=0;arm<2;arm++) {
        /* Translation-free evaluation of g*x-e. Homogeneity makes g*x zero
         * analytically; keeping the measured value avoids relying on it. */
        for(int k=1;k<3;k++)for(int axis=0;axis<2;axis++)
            intercept[arm]+=g[arm][2*k+axis]*(uv[2*(size_t)face[k]+axis]-uv[2*(size_t)face[0]+axis]);
        intercept[arm]-=residual[arm];
    }
    for(int k=0;k<3;k++) {
        size_t vertex=(size_t)face[k],edge=(size_t)system->face_edge[3*f+(size_t)k];
        int next=(k+1)%3;
        for(int arm=0;arm<2;arm++) {
            double x=g[arm][2*k],y=g[arm][2*k+1];
            double nx=g[arm][2*next],ny=g[arm][2*next+1];
            system->diagonal[3*vertex]+=d->weight*x*x;
            system->diagonal[3*vertex+1]+=d->weight*x*y;
            system->diagonal[3*vertex+2]+=d->weight*y*y;
            system->edge_block[3*edge]+=d->weight*x*nx;
            system->edge_block[3*edge+1]+=d->weight*x*ny;
            system->edge_block[3*edge+2]+=d->weight*y*ny;
            system->rhs[2*vertex]+=d->weight*x*intercept[arm];
            system->rhs[2*vertex+1]+=d->weight*y*intercept[arm];
        }
    }
    return 0;
}

static double sra_derivative(const SRDirection *d,const int32_t *face,
    const double *current,const double *target)
{
    double r[2],g[2][6],result=0;
    if(sra_model(d,face,current,r,g))return INFINITY;
    for(int arm=0;arm<2;arm++)for(int k=0;k<3;k++)for(int axis=0;axis<2;axis++) {
        size_t index=2*(size_t)face[k]+(size_t)axis;
        result+=2*d->weight*r[arm]*g[arm][2*k+axis]*(target[index]-current[index]);
    }
    return result;
}
