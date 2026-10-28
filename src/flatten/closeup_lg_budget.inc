/* Private local area-budget constraints. The exact physical acceptance rule
 * remains unchanged. When the currently good near-threshold faces in any
 * original group collectively exceed its remaining bad10 area allowance,
 * keep those faces inside the 10% bounds in this quadratic's linear model.
 * Already-bad faces may recover; freed capacity is recomputed each iteration.
 * This conservative active frontier is not an exact nonlinear area optimizer. */
typedef struct SmBudgetCurve SmBudgetCurve;
typedef struct SmBudget {
    const SmMetric *metric;size_t groups;const uint8_t *masks;double *area;
    const double *limits;
    size_t calls,rows,active_faces,active_groups,peak_faces,peak_groups;
    double seconds;int guard_all;size_t ambiguous_skipped;uint8_t *feedback;SmBudgetCurve *curvature;
} SmBudget;
static int smb_setup(Arena_T arena,const SmMetric *m,size_t groups,const uint8_t *masks,SmBudget *b)
{
    if(!arena||!m||!m->nf||!m->rest||!m->faces||!groups||!masks||!b||groups>SIZE_MAX/m->nf||groups>SIZE_MAX/sizeof(double))return -1;
    memset(b,0,sizeof *b);b->metric=m;b->groups=groups;b->masks=masks;
    b->area=ARENA_CALLOC(arena,groups,sizeof *b->area);
    for(size_t f=0;f<m->nf;f++) {
        double area=m->rest[f].area;if(!(area>0)||!isfinite(area))return -1;
        for(size_t g=0;g<groups;g++)if(masks[g*m->nf+f])b->area[g]+=area;
    }
    for(size_t g=0;g<groups;g++)if(!(b->area[g]>0)||!isfinite(b->area[g]))return -1;
    return 0;
}
static int smb_bounds(Arena_T arena,SmBudget *b,const int32_t *map,const double *uv,SmBound **out,size_t *count)
{
    if(!arena||!b||!b->metric||!b->groups||!map||!uv||!out||!count)return -1;
    double started=ves_clock_sec();const SmMetric *m=b->metric;
    if(m->nf>(SIZE_MAX-*count)/2||*count+2*m->nf>SIZE_MAX/sizeof(SmBound))return -1;
    const double lower=1/1.1,upper=1.1,skin=1e-4;
    double *sigma=ARENA_ALLOC(arena,2*m->nf*sizeof *sigma);
    double *used=ARENA_CALLOC(arena,b->groups,sizeof *used),*frontier=ARENA_CALLOC(arena,b->groups,sizeof *frontier);
    uint8_t *good_faces=ARENA_CALLOC(arena,m->nf,1);
    uint8_t *frontier_faces=ARENA_CALLOC(arena,m->nf,1),*risk=ARENA_CALLOC(arena,b->groups,1);
    for(size_t f=0;f<m->nf;f++) {
        double j[4];smm_jacobian(m->rest+f,m->faces+3*f,uv,j);
        if(smm_singular(j,sigma+2*f,NULL))return -1;
        int bad=sigma[2*f]<lower||sigma[2*f+1]>upper;
        if(m->exact_xyz){unsigned bits;if(smm_exact(m,f,uv,&bits))return -1;bad=(bits&3)!=0;}
        good_faces[f]=(uint8_t)!bad;
        frontier_faces[f]=(uint8_t)(!bad&&((b->feedback&&b->feedback[f])||sigma[2*f]<=lower+skin||sigma[2*f+1]>=upper-skin));
        for(size_t g=0;g<b->groups;g++)if(b->masks[g*m->nf+f]) {
            if(bad)used[g]+=m->rest[f].area;
            if(frontier_faces[f])frontier[g]+=m->rest[f].area;
        }
    }
    double tolerance=1e-10*fmax(1,b->area[0]);size_t active_groups=0;
    for(size_t g=0;g<b->groups;g++) {
        double capacity=fmin(.05*b->area[g],b->limits?b->limits[g]:.05*b->area[g]);
        if(used[g]>capacity+tolerance)return -1;
        risk[g]=(uint8_t)(used[g]+frontier[g]>capacity+tolerance);active_groups+=risk[g]!=0;
    }
    if(!active_groups){b->calls++;b->seconds+=ves_clock_sec()-started;return 0;}
    SmBound *rows=ARENA_ALLOC(arena,(*count+2*m->nf)*sizeof *rows);size_t n=*count,active_faces=0,skipped=0;
    if(*count)memcpy(rows,*out,*count*sizeof *rows);
    for(size_t f=0;f<m->nf;f++)if(good_faces[f]&&(b->guard_all||frontier_faces[f])) {
        int active=0;for(size_t g=0;g<b->groups;g++)if(risk[g]&&b->masks[g*m->nf+f]){active=1;break;}
        if(!active)continue;
        const int32_t *face=m->faces+3*f;double j[4],s[2],gradient[8];
        smm_jacobian(m->rest+f,face,uv,j);if(smm_singular(j,s,gradient))return -1;
        /* Far repeated singular values have no unique spectral derivative.
         * Leave these to the unchanged nonlinear guard; near-bound ambiguity
         * remains an explicit refusal, as in the original metric model. */
        if(s[1]-s[0]<1e-10*fmax(1,s[1])){if(frontier_faces[f])return -1;skipped++;continue;}
        active_faces++;
        for(int side=0;side<2;side++)if(b->guard_all||(b->feedback&&(b->feedback[f]&(1u<<side)))||(side?s[1]>=upper-skin:s[0]<=lower+skin)) {
            SmBound *r=rows+n++;memset(r,0,sizeof *r);r->face=f;r->upper=side?4:3;
            r->rhs=side?s[1]-(upper-skin):(lower+skin)-s[0];double sign=side?-1:1;
            for(size_t k=0;k<3;k++)for(size_t xy=0;xy<2;xy++) {
                if(map[face[k]]<0||map[face[k]]>INT_MAX/2-1)return -1;
                r->index[2*k+xy]=2*map[face[k]]+(int)xy;
                r->g[2*k+xy]=sign*(gradient[4*side+2*xy]*m->rest[f].gx[k]+
                                                   gradient[4*side+2*xy+1]*m->rest[f].gy[k]);
            }
        }
    }
    b->calls++;b->ambiguous_skipped+=skipped;b->rows+=n-*count;b->active_faces+=active_faces;b->active_groups+=active_groups;
    b->peak_faces=b->peak_faces>active_faces?b->peak_faces:active_faces;
    b->peak_groups=b->peak_groups>active_groups?b->peak_groups:active_groups;
    *out=rows;*count=n;b->seconds+=ves_clock_sec()-started;return 0;
}
