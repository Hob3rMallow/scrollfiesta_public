#include "closeup_local_global.h"
#include "sparse_solve.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "closeup_lg_motion.inc"
#include "closeup_lg_feedback.inc"
#include "closeup_lg_refine.inc"
#include "closeup_lg_rim_trace.inc"

typedef struct ClgContext {
    Arena_T arena;
    const float *xyz;
    const int32_t *faces;
    size_t nv,nf;
    const PMAOptions *caller;
    SmFeedback feedback;
    SmRimTrace trace;
    SmRim rim;
} ClgContext;

static int clg_material(void *context,const double *uv,size_t *pairs)
{
    ClgContext *s=context;
    return s->caller->material_guard(s->caller->context,uv,pairs);
}

static int clg_observe(void *context,const double *uv,const PMAEvent *e,
                       const PMAMetric *metric,size_t groups)
{
    ClgContext *s=context;
    if(srt_step(&s->trace,uv,e)<0||smf_observe(&s->feedback,uv,e,metric,groups))return -1;
    return s->caller->observe?s->caller->observe(s->caller->context,uv,e,metric,groups):0;
}

static int clg_refine(void *context,const double *current,const double *trial,
                      const PMAEvent *e,const PMAMetric *metric,size_t groups)
{
    (void)trial;(void)metric;(void)groups;
    ClgContext *d=context;SmFeedback *s=&d->feedback;
    Arena_Mark mark=Arena_save(d->arena);volatile int retry=-1;
    TRY {
        size_t added10=0,added25=0,rim_added=0,curve_added=0,budget_added=0;
        int pending=s->pending&&s->iteration==e->iteration&&s->refinement==e->refinement&&s->line+1==e->line;
        int rim_pending=d->trace.pending&&d->trace.previous.iteration==e->iteration&&
            d->trace.previous.refinement==e->refinement&&d->trace.previous.line==e->line-1;
        int status=0;
        if(rim_pending) {
            int32_t *witnesses=NULL,*new_pairs=NULL;size_t count=0;
            status=CloseupAngular_collect_rim(d->arena,d->xyz,d->nv,d->faces,d->nf,
                current,d->trace.refused,&witnesses,&count);
            if(!status&&count!=d->trace.previous.rim_path_failures)status=-1;
            if(!status)status=smrl_learn(d->arena,&d->rim,current,witnesses,count,&new_pairs,&rim_added);
        }
        if(!status&&pending)status=smc_learn(s->metric,current,s->pending_uv,s->pending25,&curve_added);
        if(!status&&pending)status=smbc_learn(s->budget,current,s->pending_uv,s->pending10,&budget_added);
        if(!status)retry=smf_refine(s,current,e,&added10,&added25);
        if(retry>=0&&(curve_added||budget_added||rim_added))retry=1;
    } FINALLY { Arena_restore(d->arena,mark); } END_TRY;
    return retry;
}

int CloseupLocalGlobal_run(Arena_T arena,const float *xyz,size_t nv,
    const int32_t *faces,size_t nf,const double *seed,const PMAOptions *options,
    double *output,PMAStats *stats)
{
    if(!arena||!xyz||!faces||!nv||!nf||!seed||!options||!output||!stats||
       !options->reference_uv||!options->material_guard||!options->groups||
       !options->group_masks||!options->group_count_guard||options->refinement_limit<=0)return -1;
    Arena_Mark mark=Arena_save(arena);volatile int status=-1;
    /* Independent lifetime arenas are released even if a borrowed arena throws. */
    Arena_T volatile warm_arena=NULL,learn_arena=NULL;
    memset(stats,0,sizeof *stats);memcpy(output,seed,2*nv*sizeof *output);
    TRY {
        ClgContext d={0};SmLinear motion={0};SmMetric metric={0};SmBudget budget={0};
        SmCurvature curve={0};SmBudgetCurve budget_curve={0};SmRimLearn learned={0};SmLocalWarm warm={0};
        d.arena=arena;d.xyz=xyz;d.faces=faces;d.nv=nv;d.nf=nf;d.caller=options;
        warm_arena=Arena_new_sized(1024*1024);learn_arena=Arena_new_sized(65536);
        warm.arena=warm_arena;learned.arena=learn_arena;
        if(!sm_setup(arena,xyz,nv,NULL,0,&motion)&&
           !smm_setup(arena,xyz,nv,faces,nf,options->reference_uv,&metric)&&
           !smr_setup(arena,nv,faces,nf,&d.rim)&&
           !smb_setup(arena,&metric,options->groups,options->group_masks,&budget)&&
           !smf_setup(arena,nv,&metric,&budget,options->group_count_guard,seed,&d.feedback)&&
           !smc_setup(arena,nv,&metric,seed,(size_t)options->refinement_limit,&curve)&&
           !smbc_setup(arena,nv,&budget,seed,(size_t)options->refinement_limit,&budget_curve)&&
           !srt_setup(arena,nv,seed,&d.trace)) {
            int valid=1;metric.exact_xyz=xyz;
            budget.limits=options->group_bad10_limits;
            if(budget.limits)for(size_t g=0;g<options->groups;g++)
                if(!isfinite(budget.limits[g])||budget.limits[g]<0)valid=0;
            for(size_t f=0;f<nf;f++) {
                unsigned bits=0;
                if(smm_exact(&metric,f,options->reference_uv,&bits)){valid=0;break;}
                metric.good[f]=(uint8_t)!(bits&12);
            }
            if(valid) {
                motion.metric=&metric;motion.rim=&d.rim;motion.budget=&budget;
                metric.ipm=&warm;metric.equilibrate=1;metric.origin_policy=1;
                d.rim.learning=&learned;d.feedback.accumulate=1;budget_curve.use_upper=1;
                PMAOptions o=*options;o.exact_metric_xyz=xyz;o.model=PMA_JOINT_ANGLE_V1;
                o.slim.linear_solve=sm_linear;o.slim.linear_context=&motion;
                o.context=&d;o.material_guard=clg_material;o.observe=clg_observe;o.refine=clg_refine;
                status=CloseupAngular_run(arena,xyz,nv,faces,nf,seed,&o,output,stats);
                fprintf(stderr,"[closeup local global] status=%d termination=%d accepted=%d calls=%zu warm_matches=%zu linear_seconds=%.9g residual=%.9g\n",
                    (int)status,stats->termination,stats->accepted_steps,motion.calls,warm.matched_rows,motion.seconds,motion.max_relative_residual);
            }
        }
    } FINALLY {
        Arena_T a=warm_arena,b=learn_arena;
        if(a)Arena_dispose(&a);if(b)Arena_dispose(&b);Arena_restore(arena,mark);
    } END_TRY;
    if(status){stats->status=-1;stats->termination=PMA_STOP_SOLVER_FAILURE;}
    return status;
}

#include "closeup_lg_curve_test.inc"
static void smb_test_map(double *uv,int face,double scale,double angle)
{
    double c=cos(angle),s=sin(angle),x=3*face;
    const double raw[6]={x,0,x+scale,0,x,1};
    for(int k=0;k<3;k++){uv[6*face+2*k]=c*raw[2*k]-s*raw[2*k+1];uv[6*face+2*k+1]=s*raw[2*k]+c*raw[2*k+1];}
}
#include "closeup_lg_budget_curve_test.inc"
#include "closeup_lg_budget_upper_test.inc"
#include "closeup_lg_feedback_test.inc"
#include "closeup_local_global_test.inc"

int CloseupLocalGlobal_selftest(void)
{
    return clg_test()||smc_selftest()||smbc_selftest()||smbu_selftest()||smf_selftest()||
        Sparse_local_global_selftest()||smiw_selftest()||CloseupAngular_selftest()||
        CloseupAngularRefineStep_selftest()||CloseupAngularLossStop_selftest()||
        CloseupAngularRefineBudget_selftest()||CloseupAngularGlobalRefine_selftest()||
        CloseupAngularProgressRefine_selftest();
}
