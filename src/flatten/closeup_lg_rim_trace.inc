/* Diagnostic observer only. Retain the actual committed current map and the
 * immediately preceding rim refusal for each deferred candidate. Neither a
 * deferred candidate nor an earlier refinement advances the current map. */
typedef struct SmRimTrace {
    size_t nv,exports;
    double *current,*refused;
    PMAEvent previous;
    int pending;
    FILE *ledger;
} SmRimTrace;

static int srt_setup(Arena_T arena,size_t nv,const double *seed,SmRimTrace *s)
{
    if(!arena||!nv||nv>SIZE_MAX/2/sizeof(double)||!seed||!s)return -1;
    for(size_t k=0;k<2*nv;k++)if(!isfinite(seed[k]))return -1;
    memset(s,0,sizeof *s);s->nv=nv;
    s->current=ARENA_ALLOC(arena,2*nv*sizeof *s->current);
    s->refused=ARENA_ALLOC(arena,2*nv*sizeof *s->refused);
    memcpy(s->current,seed,2*nv*sizeof *s->current);return 0;
}

/* 1 means export the captured current/refused pair; 0 means no witness.
 * The positive return is diagnostic, never an acceptance decision. */
static int srt_step(SmRimTrace *s,const double *uv,const PMAEvent *e)
{
    if(!s||!s->nv||!s->current||!s->refused||!uv||!e||e->iteration<1||e->refinement<0||e->line<0||
       (e->accepted&&e->deferred))return -1;
    int ready=e->deferred&&s->pending&&s->previous.iteration==e->iteration&&
        s->previous.refinement==e->refinement&&s->previous.line==e->line-1;
    if(e->accepted) {
        for(size_t k=0;k<2*s->nv;k++)if(!isfinite(uv[k]))return -1;
        memcpy(s->current,uv,2*s->nv*sizeof *s->current);s->pending=0;return 0;
    }
    if(e->deferred){s->pending=0;return ready;}
    s->pending=0;
    if(e->first_failure==PMA_RIM_PATH) {
        for(size_t k=0;k<2*s->nv;k++)if(!isfinite(uv[k]))return -1;
        memcpy(s->refused,uv,2*s->nv*sizeof *s->refused);
        s->previous=*e;s->pending=1;
    }
    return 0;
}

