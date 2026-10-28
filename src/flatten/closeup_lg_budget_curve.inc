/* Private directional cuts for actual active 10% area-budget rows. Original
 * group capacity, current goodness and recovery are owned by smb_bounds.
 * Upper cuts support the convex spectral norm at a measured refused trial;
 * lower cuts retain the fixed-current polar cone. Finite cuts are a model,
 * never a substitute for exact nonlinear/path/material acceptance. */
struct SmBudgetCurve {
    SmCurvature state;size_t nf,active_faces,active_peak;uint8_t *active;int use_upper;
};

static int smbc_setup(Arena_T arena,size_t nv,SmBudget *b,const double *current,
                      size_t limit,SmBudgetCurve *c)
{
    if(!arena||!b||!b->metric||!c)return -1;
    memset(c,0,sizeof *c);SmMetric local=*b->metric;
    if(smc_setup(arena,nv,&local,current,limit,&c->state))return -1;
    c->nf=local.nf;c->active=ARENA_CALLOC(arena,c->nf,1);b->curvature=c;return 0;
}

static int smbc_begin(SmBudgetCurve *c,const double *current,const SmBound *rows,size_t count)
{
    if(!c||!c->active||!current||(count&&!rows))return -1;
    for(size_t r=0;r<count;r++)if((rows[r].upper==3||(c->use_upper&&rows[r].upper==4))&&rows[r].face>=c->nf)return -1;
    if(smc_begin(&c->state,current))return -1;
    memset(c->active,0,c->nf);size_t active=0;
    for(size_t r=0;r<count;r++) {
        if(rows[r].upper==3)c->active[rows[r].face]|=1;
        if(c->use_upper&&rows[r].upper==4)c->active[rows[r].face]|=2;
    }
    for(size_t f=0;f<c->nf;f++)active+=c->active[f]!=0;
    c->active_faces+=active;c->active_peak=c->active_peak>active?c->active_peak:active;return 0;
}

static int smbc_candidate(const SmBudget *b,const double *current,const double *trial,
                           size_t face,int side,SmCurvePlane *out)
{
    const SmMetric *m=b->metric;const SmBudgetCurve *c=b->curvature;
    if(!(c->active[face]&(1u<<side)))return 0;
    double jc[4],jt[4],sc[2],st[2],gradient[8];const double lower=1/1.1;
    smm_jacobian(m->rest+face,m->faces+3*face,current,jc);
    smm_jacobian(m->rest+face,m->faces+3*face,trial,jt);
    if(smm_singular(jc,sc,gradient)||smm_singular(jt,st,NULL))return -1;
    if(m->exact_xyz){unsigned bits;if(smm_exact(m,face,current,&bits)||(bits&3))return -1;}
    else if(sc[0]<lower||sc[1]>1.1)return -1;
    double predicted=smc_dot4(gradient+4*side,jt);
    double linear_margin=side?1.1-predicted:predicted-lower;
    double actual_margin=side?1.1-st[1]:st[0]-lower;
    if(linear_margin<=1e-9||actual_margin>=-1e-9)return 0;
    SmCurvePlane candidate={0};candidate.face=face;candidate.side=side;
    if(smc_direction(jc,jt,side,candidate.gradient))return -1;
    double measured=smc_dot4(candidate.gradient,jt);
    if(side?measured<=1.1:measured>=lower)return -1;
    for(size_t i=0;i<c->state.count;i++)if(c->state.planes[i].face==face&&c->state.planes[i].side==side) {
        double error=0;for(int k=0;k<4;k++)error=fmax(error,fabs(c->state.planes[i].gradient[k]-candidate.gradient[k]));
        if(error<=1e-12)return 0;
    }
    *out=candidate;return 1;
}

static int smbc_learn(SmBudget *b,const double *current,const double *trial,
                       const uint8_t *pending,size_t *added)
{
    if(!b||!b->metric||!b->curvature||!current||!trial||!pending||!added)return -1;
    SmBudgetCurve *c=b->curvature;SmCurvature *s=&c->state;*added=0;
    if(memcmp(s->current,current,2*s->nv*sizeof *current))return -1;
    for(size_t k=0;k<2*s->nv;k++)if(!isfinite(trial[k]))return -1;
    double started=ves_clock_sec();size_t old=s->count,needed=0,skipped=0;
    for(size_t f=0;f<c->nf;f++)for(int side=0;side<(c->use_upper?2:1);side++)if(pending[f]&(1u<<side)) {
        SmCurvePlane row={0};int found=smbc_candidate(b,current,trial,f,side,&row);
        if(found<0)return -1;needed+=found!=0;skipped+=found==0;
    }
    if(s->count>s->capacity||needed>s->capacity-s->count)return -1;
    size_t write=old;
    for(size_t f=0;f<c->nf;f++)for(int side=0;side<(c->use_upper?2:1);side++)if(pending[f]&(1u<<side)) {
        SmCurvePlane row={0};int found=smbc_candidate(b,current,trial,f,side,&row);
        if(found<0)return -1;if(found)s->planes[write++]=row;
    }
    if(write!=old+needed)return -1;
    s->count=write;s->learned+=needed;s->skipped+=skipped;
    s->peak=s->peak>s->count?s->peak:s->count;s->seconds+=ves_clock_sec()-started;
    *added=needed;return 0;
}

static int smbc_bounds(Arena_T arena,SmBudget *b,const int32_t *map,
                        const double *current,SmBound **out,size_t *count)
{
    SmBudgetCurve *c=b->curvature;SmCurvature *s=c?&c->state:NULL;
    if(!s||!s->count)return 0;
    if(!arena||!map||!current||!out||!count||memcmp(s->current,current,2*s->nv*sizeof *current)||
       *count>SIZE_MAX/sizeof(SmBound)-s->count)return -1;
    const SmMetric *m=b->metric;size_t active=0;
    for(size_t i=0;i<s->count;i++)active+=(c->active[s->planes[i].face]&(1u<<s->planes[i].side))!=0;
    if(!active)return 0;
    SmBound *rows=ARENA_ALLOC(arena,(*count+active)*sizeof *rows);size_t n=*count;
    if(*count)memcpy(rows,*out,*count*sizeof *rows);
    for(size_t i=0;i<s->count;i++)if(c->active[s->planes[i].face]&(1u<<s->planes[i].side)) {
        const SmCurvePlane *p=s->planes+i;const int32_t *face=m->faces+3*p->face;
        SmBound *r=rows+n++;double j[4],sign=p->side?-1:1;
        smm_jacobian(m->rest+p->face,face,current,j);
        r->face=p->face;r->upper=p->side?4:3;
        r->rhs=sign*((p->side?1.1-1e-4:1/1.1+1e-4)-smc_dot4(p->gradient,j));
        for(size_t k=0;k<3;k++)for(size_t xy=0;xy<2;xy++) {
            r->index[2*k+xy]=2*map[face[k]]+(int)xy;
            r->g[2*k+xy]=sign*(p->gradient[2*xy]*m->rest[p->face].gx[k]+
                p->gradient[2*xy+1]*m->rest[p->face].gy[k]);
        }
    }
    *out=rows;*count=n;s->applied+=active;return 0;
}
