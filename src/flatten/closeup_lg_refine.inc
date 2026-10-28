/* Learn only newly observed spectral face/side rows from the physical failure
 * directly before this feasible candidate. Keeping current fixed allows the
 * next call to resolve the same quadratic. The nonlinear evaluator owns the
 * candidate and the final selected output throughout. */
static int smf_refine(SmFeedback *s,const double *current,const PMAEvent *e,size_t *added10,size_t *added25)
{
    *added10=*added25=0;
    if(!s||!current||!e||!s->accumulate||memcmp(current,s->current,2*s->nv*sizeof *current))return -1;
    if(!s->pending||s->iteration!=e->iteration||s->refinement!=e->refinement||s->line+1!=e->line)return 0;
    for(size_t f=0;f<s->metric->nf;f++) {
        unsigned new10=(unsigned)(s->pending10[f]&~s->budget->feedback[f]);
        unsigned new25=(unsigned)(s->pending25[f]&~s->metric->feedback[f]);
        *added10+=(new10&1)!=0;*added10+=(new10&2)!=0;
        *added25+=(new25&1)!=0;*added25+=(new25&2)!=0;
    }
    if(!*added10&&!*added25)return 0;
    for(size_t f=0;f<s->metric->nf;f++) {
        s->budget->feedback[f]|=s->pending10[f];s->metric->feedback[f]|=s->pending25[f];
    }
    return 1;
}
