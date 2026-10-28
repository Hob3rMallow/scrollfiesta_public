/* Public adapter: solver coordinates are centered and positive handed;
 * every material query uses the exact restored packed coordinates. */
#include "closeup_local_global.h"
#include "metric_bounds.h"

typedef struct CaLocalGlobal {
    CaDomain *domain;
    const CloseupAlignOpts *options;
    const double *reference;
    double center[2],*restored;
    uint8_t *original_good10;
    int orientation;
    double area[2],bad25[2];
    size_t bad25_count[2];
} CaLocalGlobal;

static int calg_budget(CaLocalGlobal *s,const double *uv,int reference)
{
    CaDomain *d=s->domain;double bad10[2]={0},bad25[2]={0},area[2]={0};size_t count[2]={0};
    const double bounds[4]={1/1.1,1.1,1/1.25,1.25};
    for(size_t f=0;f<d->nf;f++) {
        float xyz[9];double target[6];MetricBoundsResult result;
        for(size_t k=0;k<3;k++) {
            size_t v=(size_t)d->faces[3*f+k];
            memcpy(xyz+3*k,d->xyz+3*v,3*sizeof(float));
            memcpy(target+2*k,uv+2*v,2*sizeof(double));
        }
        if(MetricBounds_classify(xyz,target,bounds,&result))return 0;
        if(reference)s->original_good10[f]=(uint8_t)!(result.violated&3);
        else if(s->original_good10[f]&&(result.violated&3))return 0;
        for(size_t g=0;g<2;g++)if(!g||(d->protected_faces&&d->protected_faces[f])) {
            double a=d->face[f].area;area[g]+=a;
            if(result.violated&3)bad10[g]+=a;
            if(result.violated&12){bad25[g]+=a;count[g]++;}
        }
    }
    if(reference) {
        memcpy(s->area,area,sizeof area);memcpy(s->bad25,bad25,sizeof bad25);
        memcpy(s->bad25_count,count,sizeof count);
    }
    for(size_t g=0;g<2;g++) {
        double tolerance=1e-10*fmax(1,area[g]);
        double limit=.05*area[g];
        if(s->options->local_global_bad10_limits)limit=fmin(limit,s->options->local_global_bad10_limits[g]);
        if(bad10[g]>limit+tolerance||bad25[g]>s->bad25[g]+tolerance||
           (g&&count[g]>s->bad25_count[g]))return 0;
    }
    return 1;
}

static int calg_material(void *context,const double *uv,size_t *pairs)
{
    CaLocalGlobal *s=context;CaDomain *d=s->domain;*pairs=0;
    for(size_t v=0;v<d->nv;v++) {
        s->restored[2*v]=uv[2*v]+s->center[0];
        s->restored[2*v+1]=s->orientation*uv[2*v+1]+s->center[1];
    }
    ca_recenter(d,s->restored,s->center);
    if(!calg_budget(s,s->restored,0))return 0;
    return !s->options->candidate_guard?1:
        s->options->candidate_guard(s->options->candidate_context,d->nv,s->restored);
}

static void ca_try_local_global(CaDomain *d,const CloseupAlignOpts *opts,
    const double *original,const double *initial,const double center[2],int orientation,
    double *trial,double *best,CloseupAlignStats *report)
{
    if(!opts->local_global_iterations||report->after.fiber_rms_degrees<opts->target_rms_degrees)return;
    if(report->trials>=CLOSEUP_ALIGN_MAX_TRIALS)return;
    Arena_Mark mark=Arena_save(d->arena);
    TRY {
        CaLocalGlobal context={0};context.domain=d;context.options=opts;context.reference=original;
        memcpy(context.center,center,sizeof context.center);context.orientation=orientation;
        context.restored=ARENA_ALLOC(d->arena,2*d->nv*sizeof(double));
        context.original_good10=ARENA_ALLOC(d->arena,d->nf);
        if(calg_budget(&context,original,1)&&calg_budget(&context,best,0)) {
            size_t protected_count=0;
            if(d->protected_faces)for(size_t f=0;f<d->nf;f++)protected_count+=d->protected_faces[f]!=0;
            size_t good_group=protected_count?2:1,groups=good_group+1,good_count=0;
            uint8_t *masks=ARENA_CALLOC(d->arena,groups*d->nf,1),flags[3]={0,1,1};
            double limits[3]={.05*context.area[0],.05*context.area[1],0};
            memset(masks,1,d->nf);
            if(protected_count)memcpy(masks+d->nf,d->protected_faces,d->nf);
            if(opts->local_global_bad10_limits)for(size_t g=0;g<good_group;g++)
                limits[g]=fmin(limits[g],opts->local_global_bad10_limits[g]);
            /* Public inputs need not carry the historical source/tile masks.
             * Preserve every originally good face, which preserves bad10 area
             * on every immutable subset, including those absent from the CLI. */
            for(size_t f=0;f<d->nf;f++) {
                good_count+=(masks[good_group*d->nf+f]=context.original_good10[f]);
            }
            limits[good_group]=0;
            if(!good_count)groups=good_group;
            double *seed=ARENA_ALLOC(d->arena,2*d->nv*sizeof(double));
            for(size_t v=0;v<d->nv;v++) {
                seed[2*v]=best[2*v]-center[0];seed[2*v+1]=orientation*(best[2*v+1]-center[1]);
            }
            PMAOptions o={0};o.slim=opts->slim;o.slim.directions=d->solver_directions;
            o.slim.direction_weight=16;o.slim.direction_rigid_gauge=0;o.slim.strip_pack=0;o.slim.guard_boundary=1;
            o.model=PMA_JOINT_ANGLE_V1;o.iterations=opts->local_global_iterations;o.refinement_limit=64;
            o.loss_window=100;o.loss_relative_tolerance=1e-8;o.target_rms_degrees=opts->target_rms_degrees;
            o.reference_uv=initial;o.groups=groups;o.group_masks=masks;o.group_count_guard=flags;
            o.group_bad10_limits=limits;
            o.context=&context;o.material_guard=calg_material;
            PMAStats stats={0};report->local_global_attempted=1;
            int status=CloseupLocalGlobal_run(d->arena,d->xyz,d->nv,d->faces,d->nf,seed,&o,trial,&stats);
            report->local_global_status=status;report->local_global_termination=stats.termination;
            report->local_global_loss_converged=stats.loss_converged;
            report->local_global_initial_loss=stats.initial_loss;report->local_global_final_loss=stats.final_loss;
            CloseupAlignTrial *record=report->trial+report->trials++;
            record->seed=3;record->direction_weight=16;record->solver_status=status;
            record->iterations_run=(size_t)stats.iterations;record->accepted_steps=(size_t)stats.accepted_steps;
            /* A failed later solve may leave a fully checked earlier winner.
             * Recheck its actual returned bytes; never label that convergence. */
            size_t pairs=0;
            if(stats.accepted_steps&&calg_material(&context,trial,&pairs)>0)
                (void)ca_output_candidate(d,opts,trial,center,orientation,best,report,record);
        }
    } FINALLY { Arena_restore(d->arena,mark); } END_TRY;
}
