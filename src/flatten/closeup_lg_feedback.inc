/* Learn only from the physical rejection immediately preceding acceptance.
 * This selects spectral tangent planes for the next quadratic, never accepts
 * material or changes the nonlinear area/path/material limits. The selection
 * may be replaced or retained as a learned set. Current face status and local
 * area capacity are always recomputed when building the next quadratic. */
typedef struct SmFeedback {
    size_t nv,observations,physical_trials,commits,faces10,faces25,peak10,peak25;
    SmMetric *metric;SmBudget *budget;const uint8_t *count_guard;
    double *current,*pending_uv,*reference25,seconds;size_t *reference25_count,group;
    uint8_t *pending10,*pending25;
    int pending,iteration,line,refinement,last_used,accumulate;
} SmFeedback;

static int smf_setup(Arena_T arena,size_t nv,SmMetric *m,SmBudget *b,
                     const uint8_t *count_guard,const double *seed,SmFeedback *s)
{
    if(!arena||!nv||nv>SIZE_MAX/2/sizeof(double)||!m||!m->nf||!m->good||
       !b||b->metric!=m||!count_guard||!seed||!s)return -1;
    for(size_t k=0;k<2*nv;k++)if(!isfinite(seed[k]))return -1;
    memset(s,0,sizeof *s);s->nv=nv;s->metric=m;s->budget=b;s->count_guard=count_guard;
    s->current=ARENA_ALLOC(arena,2*nv*sizeof *s->current);memcpy(s->current,seed,2*nv*sizeof *seed);
    s->pending_uv=ARENA_ALLOC(arena,2*nv*sizeof *s->pending_uv);
    s->pending10=ARENA_CALLOC(arena,m->nf,1);s->pending25=ARENA_CALLOC(arena,m->nf,1);
    m->feedback=ARENA_CALLOC(arena,m->nf,1);b->feedback=ARENA_CALLOC(arena,m->nf,1);
    s->reference25=ARENA_CALLOC(arena,b->groups,sizeof *s->reference25);
    s->reference25_count=ARENA_CALLOC(arena,b->groups,sizeof *s->reference25_count);
    for(size_t f=0;f<m->nf;f++)if(!m->good[f])for(size_t g=0;g<b->groups;g++)if(b->masks[g*m->nf+f]) {
        s->reference25[g]+=m->rest[f].area;s->reference25_count[g]++;
    }
    return 0;
}

static int smf_observe(SmFeedback *s,const double *uv,const PMAEvent *e,const PMAMetric *metric,size_t groups)
{
    if(!s||!s->metric||!s->budget||!uv||!e||!metric||groups!=s->budget->groups)return -1;
    double started=ves_clock_sec();SmMetric *m=s->metric;SmBudget *b=s->budget;s->observations++;
    if(e->accepted) {
        for(size_t k=0;k<2*s->nv;k++)if(!isfinite(uv[k]))return -1;
        int use=s->pending&&s->iteration==e->iteration&&s->refinement==e->refinement&&s->line+1==e->line;s->last_used=use;
        size_t n10=0,n25=0;
        for(size_t f=0;f<m->nf;f++) {
            b->feedback[f]=(uint8_t)((s->accumulate?b->feedback[f]:0)|(use?s->pending10[f]:0));
            m->feedback[f]=(uint8_t)((s->accumulate?m->feedback[f]:0)|(use?s->pending25[f]:0));
            n10+=b->feedback[f]!=0;n25+=m->feedback[f]!=0;
        }
        s->commits+=use;s->faces10+=n10;s->faces25+=n25;
        s->peak10=s->peak10>n10?s->peak10:n10;s->peak25=s->peak25>n25?s->peak25:n25;
        memcpy(s->current,uv,2*s->nv*sizeof *uv);s->pending=0;
    } else {
        s->pending=0;memset(s->pending10,0,m->nf);memset(s->pending25,0,m->nf);
        if(e->first_failure==PMA_PHYSICAL) {
            size_t g=e->physical_group;if(g>=groups)return -1;
            double tolerance=1e-10*fmax(1,b->area[0]);
            double limit=fmin(.05*metric[g].area,b->limits?b->limits[g]:.05*metric[g].area);
            int bad10=metric[g].bad10_area>limit+tolerance;
            int bad25=metric[g].bad25_area>s->reference25[g]+tolerance||
                (s->count_guard[g]&&metric[g].bad25_faces>s->reference25_count[g]);
            s->physical_trials++;
            if(bad10||bad25)for(size_t f=0;f<m->nf;f++)if(b->masks[g*m->nf+f]) {
                double j[4],now[2],trial[2];smm_jacobian(m->rest+f,m->faces+3*f,s->current,j);
                if(smm_singular(j,now,NULL))return -1;
                smm_jacobian(m->rest+f,m->faces+3*f,uv,j);if(smm_singular(j,trial,NULL))return -1;
                if(m->exact_xyz) {
                    unsigned nb,tb;if(smm_exact(m,f,s->current,&nb)||smm_exact(m,f,uv,&tb))return -1;
                    if(bad10&&!(nb&3))s->pending10[f]=(uint8_t)(tb&3);
                    if(bad25&&m->good[f])s->pending25[f]=(uint8_t)((tb>>2)&3);
                } else {
                    if(bad10&&now[0]>=1/1.1&&now[1]<=1.1)
                        s->pending10[f]=(uint8_t)((trial[0]<1/1.1?1:0)|(trial[1]>1.1?2:0));
                    if(bad25&&m->good[f])s->pending25[f]=(uint8_t)((trial[0]<.8?1:0)|(trial[1]>1.25?2:0));
                }
            }
            s->pending=bad10||bad25;s->iteration=e->iteration;s->line=e->line;s->refinement=e->refinement;s->group=g;
            if(s->pending)memcpy(s->pending_uv,uv,2*s->nv*sizeof *uv);
        }
    }
    s->seconds+=ves_clock_sec()-started;return 0;
}
