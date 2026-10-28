/* Derivatives of the actual bilinear map at original measured fibre nodes.
 * A grid edge has two derivatives, and a grid vertex has four. Include every
 * incident cell, splitting the observation's physical mass across them; do
 * not invent measurements when the map basis is refined. */
typedef struct FmDirectionStencil {
    int node[4], count;
    double h[4], v[4], base_h[2], base_v[2], weight;
    double dual[4];
} FmDirectionStencil;

static int fm_direction_cells(double coordinate,int nodes,int cell[2],double fraction[2])
{
    if(!isfinite(coordinate) || coordinate<0 || coordinate>nodes-1)return 0;
    double nearest=round(coordinate),tol=16*DBL_EPSILON*fmax(1,fabs(coordinate));
    int count=0;
    if(fabs(coordinate-nearest)<=tol){
        int at=(int)nearest;
        if(at>0){cell[count]=at-1;fraction[count++]=1;}
        if(at<nodes-1){cell[count]=at;fraction[count++]=0;}
    }else{cell[0]=(int)floor(coordinate);fraction[0]=coordinate-cell[0];count=1;}
    return count;
}

static int fm_direction_stencils(Arena_T arena,const FiberMapDirections *directions,
                                  const FiberMapResult *map,const double *base,
                                  FmDirectionStencil **out,size_t *count)
{
    *out=NULL;*count=0;
    if(!directions)return 0;
    if(!directions->count || !directions->sample ||
       directions->count>SIZE_MAX/(4*sizeof(FmDirectionStencil)))return -1;
    size_t needed=0;
    for(size_t i=0;i<directions->count;i++) {
        const FiberMapDirection *p=directions->sample+i;
        int cx[2],cy[2];double fx[2],fy[2];
        int nx=fm_direction_cells((p->column-map->col0)/map->cell,map->cols,cx,fx);
        int ny=fm_direction_cells((p->row-map->row0)/map->cell,map->rows,cy,fy);
        if(!nx||!ny||needed>SIZE_MAX/sizeof(FmDirectionStencil)-(size_t)(nx*ny))return -1;
        needed+=(size_t)(nx*ny);
    }
    FmDirectionStencil *all=ARENA_CALLOC(arena,needed,sizeof *all);
    for(size_t i=0;i<directions->count;i++){
        const FiberMapDirection *p=directions->sample+i;
        double determinant=p->h[0]*p->v[1]-p->h[1]*p->v[0];
        if(!(p->weight>0) || !isfinite(p->weight) || !(determinant>0) || !isfinite(determinant))return -1;
        int cx[2],cy[2];double fx[2],fy[2];
        int nx=fm_direction_cells((p->column-map->col0)/map->cell,map->cols,cx,fx);
        int ny=fm_direction_cells((p->row-map->row0)/map->cell,map->rows,cy,fy);
        if(!nx || !ny)return -1;
        for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
            FmDirectionStencil *s=all+(*count)++;
            size_t a=(size_t)cy[y]*(size_t)map->cols+(size_t)cx[x];
            size_t node[4]={a,a+1,a+(size_t)map->cols,a+(size_t)map->cols+1};
            double gx[4]={-(1-fy[y]),1-fy[y],-fy[y],fy[y]};
            double gy[4]={-(1-fx[x]),-fx[x],1-fx[x],fx[x]};
            for(int k=0;k<4;k++){
                double h=(gx[k]*p->h[0]+gy[k]*p->h[1])/map->cell;
                double v=(gx[k]*p->v[0]+gy[k]*p->v[1])/map->cell;
                if(!isfinite(h) || !isfinite(v))return -1;
                if(h==0 && v==0)continue;
                s->node[s->count]=(int)node[k];s->h[s->count]=h;s->v[s->count++]=v;
            }
            if(s->count<2)return -1;
            s->h[0]=s->v[0]=0;
            size_t origin=2*(size_t)s->node[0];
            for(int k=1;k<s->count;k++){
                s->h[0]-=s->h[k];s->v[0]-=s->v[k];
                size_t q=2*(size_t)s->node[k];
                for(int d=0;d<2;d++){
                    double change=base[q+(size_t)d]-base[origin+(size_t)d];
                    s->base_h[d]+=s->h[k]*change;s->base_v[d]+=s->v[k]*change;
                }
            }
            s->weight=p->weight/(nx*ny);
        }
    }
    *out=all;return 0;
}

/* The measured fourth-order field has four equivalent arms. Its two families
 * are an unordered pair: assigning both arms to the same output axis is not
 * an alignment. Evaluate the better JOINT assignment, independently of the
 * local proximal targets. Original-face and cell orientation are separate. */
static double fm_direction_error(const FmDirectionStencil *stencil,size_t count,const double *uv)
{
    double sum=0,mass=0;
    for(size_t i=0;i<count;i++){
        const FmDirectionStencil *s=stencil+i;double h[2]={0},v[2]={0};
        size_t origin=2*(size_t)s->node[0];
        for(int k=1;k<s->count;k++)for(int d=0;d<2;d++){
            double change=uv[2*(size_t)s->node[k]+(size_t)d]-uv[origin+(size_t)d];
            h[d]+=s->h[k]*change;v[d]+=s->v[k]*change;
        }
        if(!isfinite(h[0]) || !isfinite(h[1]) || !isfinite(v[0]) || !isfinite(v[1]) ||
           !(hypot(h[0],h[1])>0) || !(hypot(v[0],v[1])>0))return INFINITY;
        double ah=atan2(fabs(h[1]),fabs(h[0])),av=atan2(fabs(v[0]),fabs(v[1]));
        double bh=atan2(fabs(h[0]),fabs(h[1])),bv=atan2(fabs(v[1]),fabs(v[0]));
        sum+=s->weight*fmin(ah*ah+av*av,bh*bh+bv*bv);mass+=2*s->weight;
    }
    return mass>0 && isfinite(sum) && isfinite(mass)?sum/mass:INFINITY;
}

/* Exact proximal step for the minimum of four signed axis assignments.
 * j/z are [H_u,V_u,H_v,V_v]. Along-arm lengths remain soft; cross-axis
 * components retain unit pressure. Taking the minimum AFTER minimizing each
 * quadratic makes this independent of a local quarter-turn relabeling. It is
 * not a proof of global convergence of the resulting nonconvex fit. */
static void fm_direction_prox(const double j[4],double alpha,double beta,double z[4])
{
    static const double target[4][4]={{1,0,0,1},{0,-1,1,0},{-1,0,0,-1},{0,1,-1,0}};
    double best=INFINITY,best_turn=INFINITY;
    for(int turn=0;turn<4;turn++){
        double candidate[4],cost=0;
        for(int k=0;k<4;k++){
            double t=target[turn][k],weight=t!=0?alpha:1;
            candidate[k]=(beta*j[k]+weight*t)/(beta+weight);
            double residual=candidate[k]-t,step=candidate[k]-j[k];
            cost+=weight*residual*residual+beta*step*step;
        }
        /* Equal minima at the half-quarter boundary need an output-space
         * tie-break. A target index would depend on the input arm labels.
         * The polar angle of candidate*j^T is unchanged by their relabeling. */
        double rotation=atan2(candidate[2]*j[0]+candidate[3]*j[1]-candidate[0]*j[2]-candidate[1]*j[3],
                              candidate[0]*j[0]+candidate[1]*j[1]+candidate[2]*j[2]+candidate[3]*j[3]);
        double tolerance=16*DBL_EPSILON*fmax(1,fmax(cost,best));
        if(!isfinite(best) || cost<best-tolerance || (fabs(cost-best)<=tolerance && rotation<best_turn)){
            best=cost;best_turn=rotation;memcpy(z,candidate,sizeof candidate);
        }
    }
}

static void fm_direction_rhs(FmDirectionStencil *stencil,size_t count,size_t nn,
                             const double *delta,double *rhs)
{
    double alpha=FIBER_MAP_DIRECTION_LENGTH_WEIGHT,beta=FIBER_MAP_DIRECTION_PENALTY;
    for(size_t i=0;i<count;i++){
        FmDirectionStencil *s=stencil+i;
        double j[4]={s->base_h[0],s->base_v[0],s->base_h[1],s->base_v[1]},z[4];
        size_t origin=(size_t)s->node[0];
        for(int k=1;k<s->count;k++){
            size_t q=(size_t)s->node[k];
            double u=delta[q]-delta[origin],v=delta[nn+q]-delta[nn+origin];
            j[0]+=s->h[k]*u;j[1]+=s->v[k]*u;j[2]+=s->h[k]*v;j[3]+=s->v[k]*v;
        }
        for(int k=0;k<4;k++)j[k]+=s->dual[k];
        fm_direction_prox(j,alpha,beta,z);
        for(int k=0;k<4;k++){s->dual[k]=j[k]-z[k];z[k]-=s->dual[k];}
        z[0]-=s->base_h[0];z[1]-=s->base_v[0];z[2]-=s->base_h[1];z[3]-=s->base_v[1];
        for(int k=0;k<s->count;k++){
            size_t q=(size_t)s->node[k];double weight=beta*s->weight;
            rhs[q]+=weight*(s->h[k]*z[0]+s->v[k]*z[1]);
            rhs[nn+q]+=weight*(s->h[k]*z[2]+s->v[k]*z[3]);
        }
    }
}

static double fm_fit_gain(double before,double after,int directions)
{
    /* An already aligned physical baseline is legitimate. Only arithmetic
     * zero receives this treatment; a nonzero angle still needs improvement. */
    double zero=256*DBL_EPSILON*DBL_EPSILON;
    if(directions && before<=zero)return after<=zero?1:0;
    return before>0?1-after/before:0;
}
