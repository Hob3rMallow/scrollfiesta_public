/* Private cutting planes for measured 25% curvature failures. Upper planes
 * support the convex spectral norm. Lower planes support lambda_min of
 * sym(R_current^T J), where the current polar rotation remains fixed during
 * refinement. Requiring that whole symmetric matrix >= L I would imply
 * sigma_min(J) >= L. A finite set of its directional planes is only a
 * relaxation: the original exact nonlinear/path/material checks still own
 * every candidate and final acceptance. */
typedef struct SmCurvePlane {
    size_t face; int side; double gradient[4];
} SmCurvePlane;
struct SmCurvature {
    size_t nv,capacity,count,learned,peak,applied,resets,skipped;
    double *current,seconds; SmCurvePlane *planes;
};

static double smc_dot4(const double *a,const double *b)
{
    double value=0;for(int k=0;k<4;k++)value+=a[k]*b[k];return value;
}

/* Output stays unchanged on refusal. The lower direction is an eigenvector
 * of the symmetric polar-frame trial, not a lower singular-value tangent at
 * that trial (which would not describe the same fixed inner cone). */
static int smc_direction(const double *current,const double *trial,int side,double *gradient)
{
    if(!current||!trial||!gradient||side<0||side>1)return -1;
    double sc[2],st[2],g[8],out[4];
    if(smm_singular(current,sc,NULL)||smm_singular(trial,st,side?g:NULL))return -1;
    if(side)memcpy(out,g+4,sizeof out);
    else {
        double a=current[0]+current[3],b=current[2]-current[1],length=hypot(a,b);
        if(!(length>0)||!isfinite(length))return -1;
        double c=a/length,s=b/length;
        double xx=c*trial[0]+s*trial[2],yy=-s*trial[1]+c*trial[3];
        double xy=.5*(c*trial[1]+s*trial[3]-s*trial[0]+c*trial[2]);
        double angle=.5*atan2(2*xy,xx-yy),x=-sin(angle),y=cos(angle);
        double rx=c*x-s*y,ry=s*x+c*y;
        out[0]=rx*x;out[1]=rx*y;out[2]=ry*x;out[3]=ry*y;
    }
    for(int k=0;k<4;k++)if(!isfinite(out[k]))return -1;
    memcpy(gradient,out,sizeof out);return 0;
}

static int smc_setup(Arena_T arena,size_t nv,SmMetric *m,const double *current,
                     size_t refinement_limit,SmCurvature *c)
{
    if(!arena||!nv||nv>SIZE_MAX/2/sizeof(double)||!m||!m->nf||!m->rest||!m->good||
       !current||!c||!refinement_limit||refinement_limit==SIZE_MAX||
       m->nf>SIZE_MAX/2/(refinement_limit+1)/sizeof(SmCurvePlane))return -1;
    for(size_t k=0;k<2*nv;k++)if(!isfinite(current[k]))return -1;
    memset(c,0,sizeof *c);c->nv=nv;c->capacity=2*m->nf*(refinement_limit+1);
    c->planes=ARENA_ALLOC(arena,c->capacity*sizeof *c->planes);
    c->current=ARENA_ALLOC(arena,2*nv*sizeof *c->current);
    memcpy(c->current,current,2*nv*sizeof *current);m->curvature=c;return 0;
}

static int smc_begin(SmCurvature *c,const double *current)
{
    if(!c||!current||!c->current||!c->nv)return -1;
    for(size_t k=0;k<2*c->nv;k++)if(!isfinite(current[k]))return -1;
    if(memcmp(c->current,current,2*c->nv*sizeof *current)) {
        c->count=0;c->resets++;memcpy(c->current,current,2*c->nv*sizeof *current);
    }
    return 0;
}

static int smc_candidate(const SmMetric *m,const double *current,const double *trial,
                         size_t face,int side,SmCurvePlane *out)
{
    double jc[4],jt[4],sc[2],st[2],gradient[8];
    smm_jacobian(m->rest+face,m->faces+3*face,current,jc);
    smm_jacobian(m->rest+face,m->faces+3*face,trial,jt);
    if(smm_singular(jc,sc,gradient)||smm_singular(jt,st,NULL))return -1;
    int included=(m->feedback&&(m->feedback[face]&(1u<<side)))||
        (side?sc[1]>=1.25-1e-4:sc[0]<=.8+1e-4);
    double predicted=smc_dot4(gradient+4*side,jt);
    double predicted_margin=side?1.25-predicted:predicted-.8;
    double actual_margin=side?1.25-st[1]:st[0]-.8;
    if(!included||predicted_margin<=1e-9||actual_margin>=-1e-9)return 0;
    SmCurvePlane candidate={0};candidate.face=face;candidate.side=side;
    if(smc_direction(jc,jt,side,candidate.gradient))return -1;
    double measured=smc_dot4(candidate.gradient,jt);
    /* A new direction must actually separate the observed nonlinear failure.
     * The lower cone can be more restrictive than the original scale bound. */
    if(side?measured<=1.25:measured>=.8)return -1;
    const SmCurvature *c=m->curvature;
    for(size_t i=0;i<c->count;i++)if(c->planes[i].face==face&&c->planes[i].side==side) {
        double error=0;for(int k=0;k<4;k++)error=fmax(error,fabs(c->planes[i].gradient[k]-candidate.gradient[k]));
        if(error<=1e-12)return 0;
    }
    *out=candidate;return 1;
}

/* Two passes make validation/capacity failure atomic. One face/side occurs at
 * most once per call; candidate construction is deterministic and only reads
 * the old set until the second pass commits the complete batch. */
static int smc_learn(SmMetric *m,const double *current,const double *trial,
                     const uint8_t *pending,size_t *added)
{
    if(!m||!m->curvature||!current||!trial||!pending||!added)return -1;
    SmCurvature *c=m->curvature;*added=0;
    if(memcmp(c->current,current,2*c->nv*sizeof *current))return -1;
    for(size_t k=0;k<2*c->nv;k++)if(!isfinite(trial[k]))return -1;
    double started=ves_clock_sec();size_t old_count=c->count,needed=0,skipped=0;
    for(size_t f=0;f<m->nf;f++)if(m->good[f])for(int side=0;side<2;side++)if(pending[f]&(1u<<side)) {
        SmCurvePlane row={0};int found=smc_candidate(m,current,trial,f,side,&row);
        if(found<0)return -1;
        needed+=found!=0;skipped+=found==0;
    }
    if(c->count>c->capacity||needed>c->capacity-c->count)return -1;
    size_t write=old_count;
    for(size_t f=0;f<m->nf;f++)if(m->good[f])for(int side=0;side<2;side++)if(pending[f]&(1u<<side)) {
        SmCurvePlane row={0};int found=smc_candidate(m,current,trial,f,side,&row);
        if(found<0)return -1;
        if(found)c->planes[write++]=row;
    }
    if(write!=old_count+needed)return -1;
    c->count=write;c->learned+=needed;c->skipped+=skipped;
    c->peak=c->peak>c->count?c->peak:c->count;c->seconds+=ves_clock_sec()-started;
    *added=needed;return 0;
}

static int smc_bounds(Arena_T arena,SmMetric *m,const int32_t *map,
                      const double *current,SmBound **out,size_t *count)
{
    SmCurvature *c=m->curvature;
    if(!c||!c->count)return 0;
    if(!arena||!map||!current||!out||!count||memcmp(c->current,current,2*c->nv*sizeof *current)||
       *count>SIZE_MAX/sizeof(SmBound)-c->count)return -1;
    SmBound *rows=ARENA_ALLOC(arena,(*count+c->count)*sizeof *rows);
    if(*count)memcpy(rows,*out,*count*sizeof *rows);
    for(size_t i=0;i<c->count;i++) {
        const SmCurvePlane *p=c->planes+i;const int32_t *face=m->faces+3*p->face;
        SmBound *r=rows+*count+i;double j[4],sign=p->side?-1:1;
        smm_jacobian(m->rest+p->face,face,current,j);
        r->face=p->face;r->upper=p->side;
        r->rhs=sign*((p->side?1.25-1e-4:.8+1e-4)-smc_dot4(p->gradient,j));
        for(size_t k=0;k<3;k++)for(size_t xy=0;xy<2;xy++) {
            r->index[2*k+xy]=2*map[face[k]]+(int)xy;
            r->g[2*k+xy]=sign*(p->gradient[2*xy]*m->rest[p->face].gx[k]+
                p->gradient[2*xy+1]*m->rest[p->face].gy[k]);
        }
    }
    *out=rows;*count+=c->count;c->applied+=c->count;return 0;
}
