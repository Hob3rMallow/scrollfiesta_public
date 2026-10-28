/* A large accepted fraction can justify a fresh local model. Near the loss
 * tolerance, keep resolving constraints so this shortcut cannot establish
 * stationarity. The exact accepted-map gates and retained best remain final. */
static int pma_progress_restart(const PMAOptions *o,double alpha,double objective,double area,double derivative)
{
    if(!(o->refinement_progress_alpha>0)||!isfinite(alpha)||!isfinite(objective)||!isfinite(area)||!isfinite(derivative))return 0;
    return alpha>=o->refinement_progress_alpha&&fabs(derivative)>o->loss_relative_tolerance*fmax(1,objective-4*area);
}

/* Refine constraints while the same SLIM/GN quadratic and current map remain
 * fixed. Every proposed candidate passes the unchanged complete evaluator.
 * The best candidate is rechecked before one final commit; later numerical,
 * no-step or refinement failures cannot discard an earlier verified candidate.
 * The limit bounds solver work, never establishes geometric infeasibility. */
static int pma_refine_step(Arena_T arena,SRSystem *s,double *current,PMAMetric *current_metric,
                           const PMAMetric *reference,const PMAOptions *o,SREdgeBox *boxes,
                           int iteration,double *output,PMAStats *stats)
{
    Arena_Mark mark=Arena_save(arena);int result=0,have_best=0,last_attempt=0;
    double *destination=ARENA_ALLOC(arena,2*s->nv*sizeof *destination);
    double *trial=ARENA_ALLOC(arena,2*s->nv*sizeof *trial),*gradient=ARENA_ALLOC(arena,2*s->nv*sizeof *gradient);
    double *best=ARENA_ALLOC(arena,2*s->nv*sizeof *best);
    PMAMetric *metric=ARENA_ALLOC(arena,o->groups*sizeof *metric);PMAEvent best_event={0};
    double best_value=DBL_MAX;
    int limit=o->refinement_limit;
    if(iteration==o->iterations&&o->final_refinement_limit>limit)limit=o->final_refinement_limit;
    for(int attempt=0;attempt<=limit;attempt++) {
        double attempt_started=ves_clock_sec();
        last_attempt=attempt;stats->refinement_solves++;
        memcpy(destination,current,2*s->nv*sizeof *destination);
        if(o->slim.linear_solve(arena,o->slim.linear_context,s->nv,s->ne,s->diagonal,s->edge_a,s->edge_b,s->edge_block,s->rhs,destination)) {
            stats->refinement_failures++;if(!have_best)result=-1;break;
        }
        double linear_elapsed=ves_clock_sec()-attempt_started;
        double derivative=pma_proxy_derivative(s,current,destination,gradient),maxstep=0;
        for(size_t k=0;k<2*s->nv;k++)maxstep=fmax(maxstep,fabs(destination[k]-current[k]));
        double h=1e-5/fmax(1,maxstep),fd=0;int finite=isfinite(derivative);
        for(size_t k=0;k<2*s->nv;k++)trial[k]=current[k]+h*(destination[k]-current[k]);
        if(pma_measure(s,trial,o,metric))finite=0;else fd=pma_objective(metric,o->model);
        for(size_t k=0;k<2*s->nv;k++)trial[k]=current[k]-h*(destination[k]-current[k]);
        if(pma_measure(s,trial,o,metric))finite=0;else fd=(fd-pma_objective(metric,o->model))/(2*h);
        if(!finite||!isfinite(fd)){stats->refinement_failures++;if(!have_best)result=-1;break;}
        int found=0,retry=0;double alpha=1;
        for(int line=0;line<o->slim.line_search;line++,alpha*=.5) {
            for(size_t k=0;k<2*s->nv;k++)trial[k]=current[k]+alpha*(destination[k]-current[k]);
            PMAEvent e={0};e.iteration=iteration;e.line=line;e.alpha=alpha;e.refinement=attempt;
            e.previous_objective=pma_objective(current_metric,o->model);e.proxy_derivative=derivative;e.fd_derivative=fd;
            found=pma_evaluate(s,current,trial,o,reference,metric,boxes,&e);stats->trials++;stats->stop_reason=e.first_failure;
            if(found) {
                if(e.objective<best_value){memcpy(best,trial,2*s->nv*sizeof *best);best_value=e.objective;best_event=e;have_best=1;}
                /* The callback may add rows, but cannot accept geometry. */
                retry=o->refine(o->context,current,trial,&e,metric,o->groups);
                if(retry<0){stats->refinement_failures++;retry=0;}
                e.accepted=0;e.deferred=1;
            }
            if(o->observe(o->context,trial,&e,metric,o->groups)){result=-1;goto done;}
            if(found)break;
        }
        if(iteration==o->iterations&&o->final_refinement_limit>0)
            fprintf(stderr,"[refinement timing] iteration=%d attempt=%d linear_seconds=%.17g seconds=%.17g found=%d retry=%d best=%.17g\n",
                iteration,attempt,linear_elapsed,ves_clock_sec()-attempt_started,found,retry,best_value);
        if(!found){stats->refinement_failures++;break;}
        if(!retry)break;
        if(pma_progress_restart(o,best_event.alpha,best_event.previous_objective,current_metric[0].area,best_event.proxy_derivative)) {
            stats->refinement_progress_stops++;
            fprintf(stderr,"[refinement progress] iteration=%d attempt=%d alpha=%.17g best=%.17g\n",iteration,attempt,best_event.alpha,best_value);
            break;
        }
        if(attempt==limit){stats->refinement_limit_hits++;break;}
        stats->refinement_retries++;
    }
    if(have_best) {
        PMAEvent e=best_event;e.accepted=0;e.deferred=0;
        int accepted=pma_evaluate(s,current,best,o,reference,metric,boxes,&e);stats->trials++;stats->stop_reason=e.first_failure;
        if(o->observe(o->context,best,&e,metric,o->groups)){result=-1;goto done;}
        if(accepted) {
            memcpy(current,best,2*s->nv*sizeof *current);memcpy(output,current,2*s->nv*sizeof *output);
            memcpy(current_metric,metric,o->groups*sizeof *metric);stats->accepted_steps++;stats->final_objective=e.objective;stats->last_proxy_derivative=e.proxy_derivative;
            stats->refinement_prior_best+=e.refinement<last_attempt;result=1;
        }
    }
    fprintf(stderr,"[private refined step] step=%d accepted=%d attempts=%d best_refinement=%d objective=%.12g RMS=%.9g stop=%d\n",
        iteration,result>0,last_attempt+1,best_event.refinement,stats->final_objective,
        current_metric[0].mass>0?sqrt(current_metric[0].angular16/(32*current_metric[0].mass))*180/3.14159265358979323846:NAN,stats->stop_reason);
done:
    Arena_restore(arena,mark);return result;
}
