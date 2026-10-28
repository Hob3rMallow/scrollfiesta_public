static int cq_axis(const double *uv,const SREdgeBox *a,const SREdgeBox *b,double normal[2],double *gap)
{
    const double *p=uv+2*(size_t)a->a,*q=uv+2*(size_t)a->b;
    const double *r=uv+2*(size_t)b->a,*t=uv+2*(size_t)b->b;
    double edges[2][2]={{q[0]-p[0],q[1]-p[1]},{t[0]-r[0],t[1]-r[1]}};
    double best=0;
    for(int e=0;e<2;e++)for(int perpendicular=0;perpendicular<2;perpendicular++) {
        double nx=perpendicular?-edges[e][1]:edges[e][0],ny=perpendicular?edges[e][0]:edges[e][1];
        double length=hypot(nx,ny);if(!(length>0))continue;nx/=length;ny/=length;
        double aq=(q[0]-p[0])*nx+(q[1]-p[1])*ny;
        double br=(r[0]-p[0])*nx+(r[1]-p[1])*ny,bt=(t[0]-p[0])*nx+(t[1]-p[1])*ny;
        double forward=fmin(br,bt)-fmax(0,aq),reverse=fmin(0,aq)-fmax(br,bt);
        if(forward>best){best=forward;normal[0]=nx;normal[1]=ny;}
        if(reverse>best){best=reverse;normal[0]=-nx;normal[1]=-ny;}
    }
    *gap=best;return best>0;
}

/* Orientation of three linearly moving points is quadratic. Prove its sign
 * on the WHOLE interval, with an upper bound on sr_orientation's unchanged
 * relative product tolerance. This permits coherent rotations/translations
 * that a fixed initial separating line cannot represent. */
static int cq_side_interval(const double *current,const double *target,int32_t ia,int32_t ib,int32_t ic)
{
    size_t a=2*(size_t)ia,b=2*(size_t)ib,c=2*(size_t)ic;
    double e[2],g[2],de[2],dg[2],em[2],gm[2];
    for(int k=0;k<2;k++) {
        e[k]=current[b+k]-current[a+k];g[k]=current[c+k]-current[a+k];
        double end_e=target[b+k]-target[a+k],end_g=target[c+k]-target[a+k];
        de[k]=end_e-e[k];dg[k]=end_g-g[k];em[k]=fmax(fabs(e[k]),fabs(end_e));gm[k]=fmax(fabs(g[k]),fabs(end_g));
    }
    double qa=de[0]*dg[1]-de[1]*dg[0];
    double qb=de[0]*g[1]-de[1]*g[0]+e[0]*dg[1]-e[1]*dg[0];
    double qc=e[0]*g[1]-e[1]*g[0],end=qa+qb+qc;
    double minimum=fmin(qc,end),maximum=fmax(qc,end);
    if(qa!=0) {double t=-.5*qb/qa;if(t>0 && t<1){double value=(qa*t+qb)*t+qc;minimum=fmin(minimum,value);maximum=fmax(maximum,value);}}
    double epsilon=1e-12*(em[0]*gm[1]+em[1]*gm[0]+1);
    if(!isfinite(epsilon) || !isfinite(minimum) || !isfinite(maximum))return 0;
    return minimum>epsilon?1:maximum<-epsilon?-1:0;
}

static int cq_swept_pair_safe(const double *current,const double *target,const SREdgeBox *a,const SREdgeBox *b)
{
    int first=cq_side_interval(current,target,a->a,a->b,b->a);
    if(first && first==cq_side_interval(current,target,a->a,a->b,b->b))return 1;
    first=cq_side_interval(current,target,b->a,b->b,a->a);
    return first && first==cq_side_interval(current,target,b->a,b->b,a->b);
}

