/* Private whole-patch uniform SLIM/GN experiment. The only linear solver is
 * the existing native TAUCS callback. A positive GN matrix is a touching local
 * model, NOT a majorizer; exact nonlinear line search owns acceptance. */
#include "closeup_angular_metric_angular.h"
#include "metric_bounds.h"
#include "closeup_angular_joint_angle.inc"
#include "closeup_angular_joint_angle_test.inc"
static const Except_T pma_failed={"private metric/angular diagnostic failure"};

static int pma_measure(SRSystem *s,const double *uv,const PMAOptions *o,PMAMetric *out)
{
    memset(out,0,o->groups*sizeof *out);
    for(size_t g=0;g<o->groups;g++)out[g].sigma_min=DBL_MAX;
    for(size_t f=0;f<s->nf;f++) {
        double j[4],angle=0;sr_face_jacobian(s->face+f,s->faces+3*f,uv,j);
        double det=j[0]*j[3]-j[1]*j[2],norm=0;
        for(int k=0;k<4;k++)norm+=j[k]*j[k];
        if(!isfinite(det)||!isfinite(norm)||!(fabs(det)>0))return -1;
        double a=j[0]*j[0]+j[1]*j[1],b=j[0]*j[2]+j[1]*j[3],c=j[2]*j[2]+j[3]*j[3];
        double large=sqrt(fmax(0,.5*(a+c+hypot(a-c,2*b))));
        double small=large>0?fabs(det)/large:0;
        if(!(small>0)||!isfinite(large))return -1;
        int bad10=small<1/1.1||large>1.1,bad25=small<1/1.25||large>1.25;
        if(o->exact_metric_xyz) {
            float physical[9];double mapped_uv[6];MetricBoundsResult bounds;
            for(int corner=0;corner<3;corner++) {
                size_t vertex=(size_t)s->faces[3*f+corner];
                memcpy(physical+3*corner,o->exact_metric_xyz+3*vertex,3*sizeof(float));
                memcpy(mapped_uv+2*corner,uv+2*vertex,2*sizeof(double));
            }
            const double limits[4]={1/1.1,1.1,1/1.25,1.25};
            if(MetricBounds_classify(physical,mapped_uv,limits,&bounds))return -1;
            bad10=(bounds.violated&3)!=0;bad25=(bounds.violated&12)!=0;
        }
        double area=s->face[f].area,sd=area*norm*(1+1/(det*det));
        const SRDirection *d=s->direction+f;
        double mapped=sr_direction_energy(d,s->faces+3*f,uv,&angle);
        double angular=sra_energy(d,s->faces+3*f,uv);
        if(!isfinite(sd)||!isfinite(mapped)||!isfinite(angular))return -1;
        for(size_t g=0;g<o->groups;g++)if(o->group_masks[g*s->nf+f]) {
            PMAMetric *m=out+g;m->faces++;m->area+=area;m->sd+=sd;
            m->mapped16+=mapped;m->angular16+=angular;m->flips+=det<=0;
            m->sigma_min=fmin(m->sigma_min,small);m->sigma_max=fmax(m->sigma_max,large);
            m->scale_error+=area*.5*((small-1)*(small-1)+(large-1)*(large-1));
            if(d->weight>0){m->measured_faces++;m->measured_area+=area;m->mass+=d->weight/16;}
            if(bad10){m->bad10_faces++;m->bad10_area+=area;}
            if(bad25){m->bad25_faces++;m->bad25_area+=area;}
        }
    }
    return 0;
}
static double pma_objective(const PMAMetric *m,int model)
{return m->sd+(model==PMA_JOINT_ANGLE_V1?m->angular16:m->mapped16);}
static int pma_physical(const PMAMetric *m,const PMAMetric *reference,const PMAOptions *o,size_t *failed)
{
    double tolerance=1e-10*fmax(1,reference[0].area);
    for(size_t g=0;g<o->groups;g++) {
        if((o->group_bad10_limits&&(!isfinite(o->group_bad10_limits[g])||o->group_bad10_limits[g]<0)) ||
           !m[g].faces || m[g].area!=reference[g].area || m[g].flips || !isfinite(m[g].sd) ||
           m[g].bad10_area>fmin(.05*m[g].area,o->group_bad10_limits?o->group_bad10_limits[g]:.05*m[g].area)+tolerance ||
           m[g].bad25_area>reference[g].bad25_area+tolerance ||
           (o->group_count_guard[g] && m[g].bad25_faces>reference[g].bad25_faces)) {
            *failed=g;return 0;
        }
    }
    *failed=SIZE_MAX;return 1;
}

static int pma_assemble(SRSystem *s,const double *uv,const PMAOptions *o)
{
    /* Both models use the identical pure-SD translation regularizer. Its
     * derivative at current is zero; it pins no material vertex. */
    SRDirection *directions=s->direction;s->direction=NULL;
    int status=sr_assemble(s,uv,&o->slim);s->direction=directions;if(status)return -1;
    for(size_t f=0;f<s->nf;f++)if(directions[f].weight>0) {
        if(o->model==PMA_JOINT_ANGLE_V1){if(sra_assemble_face(s,f,uv))return -1;}
        else {
            const SRDirection *d=directions+f;const int32_t *face=s->faces+3*f;
            int swapped=sr_direction_assignment(d,face,uv);
            const double *u=swapped?d->h:d->v,*v=swapped?d->v:d->h;
            for(size_t k=0;k<3;k++) {
                size_t vertex=(size_t)face[k],next=(k+1)%3,edge=(size_t)s->face_edge[3*f+k];
                s->diagonal[3*vertex]+=d->weight*u[k]*u[k];s->diagonal[3*vertex+2]+=d->weight*v[k]*v[k];
                s->edge_block[3*edge]+=d->weight*u[k]*u[next];s->edge_block[3*edge+2]+=d->weight*v[k]*v[next];
            }
        }
    }
    return 0;
}
static double pma_proxy_derivative(const SRSystem *s,const double *uv,const double *target,double *gradient)
{
    for(size_t v=0;v<s->nv;v++) {
        gradient[2*v]=s->diagonal[3*v]*uv[2*v]+s->diagonal[3*v+1]*uv[2*v+1]-s->rhs[2*v];
        gradient[2*v+1]=s->diagonal[3*v+1]*uv[2*v]+s->diagonal[3*v+2]*uv[2*v+1]-s->rhs[2*v+1];
    }
    for(size_t e=0;e<s->ne;e++) {
        size_t a=(size_t)s->edge_a[e],b=(size_t)s->edge_b[e];const double *w=s->edge_block+3*e;
        gradient[2*a]+=w[0]*uv[2*b]+w[1]*uv[2*b+1];gradient[2*a+1]+=w[1]*uv[2*b]+w[2]*uv[2*b+1];
        gradient[2*b]+=w[0]*uv[2*a]+w[1]*uv[2*a+1];gradient[2*b+1]+=w[1]*uv[2*a]+w[2]*uv[2*a+1];
    }
    double result=0;for(size_t k=0;k<2*s->nv;k++)result+=2*gradient[k]*(target[k]-uv[k]);return result;
}

/* Reuse the previously tested quadratic moving-edge orientation certificate.
 * Relative separating axes also prove collinear-disjoint intervals. This is
 * geometric path safety; native endpoint tolerances remain an extra gate. */
#include "closeup_angular_swept_rim.inc"
static size_t pma_rim_path_trace(const SRSystem *s,const double *uv,const double *target,SREdgeBox *boxes,FILE *trace,const int32_t *aliases,int32_t *pairs,size_t capacity)
{
    for(size_t k=0;k<s->nb;k++) {
        const SRBoundaryEdge *e=s->boundary+k;SREdgeBox *b=boxes+k;size_t a=(size_t)e->a,c=(size_t)e->b;
        b->a=e->a;b->b=e->b;b->component=0; /* all components participate */
        b->min_u=fmin(fmin(uv[2*a],uv[2*c]),fmin(target[2*a],target[2*c]));
        b->max_u=fmax(fmax(uv[2*a],uv[2*c]),fmax(target[2*a],target[2*c]));
        b->min_v=fmin(fmin(uv[2*a+1],uv[2*c+1]),fmin(target[2*a+1],target[2*c+1]));
        b->max_v=fmax(fmax(uv[2*a+1],uv[2*c+1]),fmax(target[2*a+1],target[2*c+1]));
    }
    qsort(boxes,s->nb,sizeof *boxes,sr_compare_box);size_t unproved=0;
    for(size_t i=0;i<s->nb;i++) {
        const SREdgeBox *a=boxes+i;double epsilon=1e-10*(hypot(a->max_u-a->min_u,a->max_v-a->min_v)+1);
        for(size_t k=i+1;k<s->nb;k++) {
            const SREdgeBox *b=boxes+k;if(b->min_u>a->max_u+epsilon)break;
            if(b->max_v<a->min_v-epsilon || b->min_v>a->max_v+epsilon)continue;
            if(a->a==b->a||a->a==b->b||a->b==b->a||a->b==b->b)continue;
            if(cq_swept_pair_safe(uv,target,a,b))continue;
            double normal[2],gap;int proven=0;
            if(cq_axis(uv,a,b,normal,&gap)) {
                double minimum=DBL_MAX;int32_t aa[2]={a->a,a->b},bb[2]={b->a,b->b};
                for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
                    size_t p=(size_t)aa[u],q=(size_t)bb[v];
                    minimum=fmin(minimum,(target[2*q]-target[2*p])*normal[0]+(target[2*q+1]-target[2*p+1])*normal[1]);
                }
                proven=minimum>0;
            }
            if(!proven) {
                if(pairs&&unproved<capacity) {
                    int32_t vertices[4]={a->a,a->b,b->a,b->b};memcpy(pairs+4*unproved,vertices,sizeof vertices);
                }
                unproved++;
                if(trace) {
                    int32_t vertices[4]={a->a,a->b,b->a,b->b};double normal[2],gap=NAN;
                    int axis=cq_axis(uv,a,b,normal,&gap);
                    int current=sr_segments_intersect(uv+2*(size_t)a->a,uv+2*(size_t)a->b,uv+2*(size_t)b->a,uv+2*(size_t)b->b);
                    int endpoint=sr_segments_intersect(target+2*(size_t)a->a,target+2*(size_t)a->b,target+2*(size_t)b->a,target+2*(size_t)b->b);
                    for(int k=0;k<4;k++)fprintf(trace,"%d,",vertices[k]);
                    for(int k=0;k<4;k++)fprintf(trace,"%d,",aliases?aliases[vertices[k]]:vertices[k]);
                    fprintf(trace,"%d,%d,%d,%d,%.17g,%d",s->vertex_component[a->a],s->vertex_component[b->a],current,endpoint,gap,axis);
                    for(int step=0;step<2;step++)for(int k=0;k<4;k++)for(int xy=0;xy<2;xy++)
                        fprintf(trace,",%.17g",(step?target:uv)[2*(size_t)vertices[k]+(size_t)xy]);
                    fputc('\n',trace);
                }
            }
        }
    }
    return unproved;
}
static size_t pma_rim_path(const SRSystem *s,const double *uv,const double *target,SREdgeBox *boxes)
{return pma_rim_path_trace(s,uv,target,boxes,NULL,NULL,NULL,0);}

int CloseupAngular_trace(Arena_T arena,const float *xyz,size_t nv,const int32_t *faces,size_t nf,
    const double *current,const double *target,const int32_t *aliases,const char *path)
{
    Arena_Mark mark=Arena_save(arena);SlimRefineOpts opts;SlimRefine_defaults(&opts);SRSystem system;
    if(sr_build_system(arena,xyz,nv,faces,nf,&opts,1,&system)){Arena_restore(arena,mark);return -1;}
    SREdgeBox *boxes=ARENA_ALLOC(arena,(system.nb?system.nb:1)*sizeof *boxes);
    FILE *fp=fopen(path,"wb");if(!fp){Arena_restore(arena,mark);return -1;}
    fputs("a0,a1,b0,b1,alias_a0,alias_a1,alias_b0,alias_b1,component_a,component_b,current_intersect,target_intersect,current_gap,axis_defined",fp);
    for(int step=0;step<2;step++)for(int k=0;k<4;k++)for(int xy=0;xy<2;xy++)fprintf(fp,",s%d_v%d_%c",step,k,xy?'v':'u');
    fputc('\n',fp);
    size_t unproved=pma_rim_path_trace(&system,current,target,boxes,fp,aliases,NULL,0);
    int rc=ferror(fp);rc|=fclose(fp)!=0;
    fprintf(stderr,"[source rim trace] unproved=%zu %s\n",unproved,path);
    Arena_restore(arena,mark);return rc?-1:0;
}

int CloseupAngular_collect_rim(Arena_T arena,const float *xyz,size_t nv,const int32_t *faces,size_t nf,
    const double *current,const double *target,int32_t **pairs,size_t *count)
{
    if(!arena||!xyz||!nv||!faces||!nf||!current||!target||!pairs||!count)return -1;
    *pairs=NULL;*count=0;SlimRefineOpts opts;SlimRefine_defaults(&opts);SRSystem system;
    if(sr_build_system(arena,xyz,nv,faces,nf,&opts,1,&system))return -1;
    SREdgeBox *boxes=ARENA_ALLOC(arena,(system.nb?system.nb:1)*sizeof *boxes);
    size_t n=pma_rim_path_trace(&system,current,target,boxes,NULL,NULL,NULL,0);
    if(n>SIZE_MAX/4/sizeof(int32_t))return -1;
    if(!n)return 0;
    int32_t *out=ARENA_ALLOC(arena,4*n*sizeof *out);
    if(pma_rim_path_trace(&system,current,target,boxes,NULL,NULL,out,n)!=n)return -1;
    *pairs=out;*count=n;return 0;
}

static int pma_evaluate(SRSystem *s,const double *current,const double *trial,const PMAOptions *o,
    const PMAMetric *reference,PMAMetric *metric,SREdgeBox *boxes,PMAEvent *e)
{
    e->first_failure=PMA_OK;e->physical_group=SIZE_MAX;
    if(pma_measure(s,trial,o,metric)){e->first_failure=PMA_NONFINITE;return 0;}
    e->objective=pma_objective(metric,o->model);
    for(size_t f=0;f<s->nf;f++)if(!UvGuard_interval(current,trial,s->faces+3*f,1))e->face_path_failures++;
    if(e->face_path_failures){e->first_failure=PMA_FACE_PATH;return 0;}
    e->rim_path_failures=pma_rim_path(s,current,trial,boxes);
    if(e->rim_path_failures){e->first_failure=PMA_RIM_PATH;return 0;}
    e->native_boundary_contacts=sr_boundary_intersections(s,trial,boxes);
    if(e->native_boundary_contacts){e->first_failure=PMA_NATIVE_BOUNDARY;return 0;}
    if(!(e->objective<e->previous_objective)){e->first_failure=PMA_OBJECTIVE;return 0;}
    if(!pma_physical(metric,reference,o,&e->physical_group)){e->first_failure=PMA_PHYSICAL;return 0;}
    int material=o->material_guard(o->context,trial,&e->material_pairs);
    if(material<=0){e->first_failure=material<0?PMA_MATERIAL_UNKNOWN:PMA_MATERIAL;return 0;}
    e->accepted=1;return 1;
}

#include "closeup_angular_refine_step.inc"
#include "closeup_angular_loss_stop.inc"

int CloseupAngular_run(Arena_T arena,const float *xyz,size_t nv,const int32_t *faces,size_t nf,
    const double *seed,const PMAOptions *o,double *output,PMAStats *stats)
{
    if(!arena||!xyz||!faces||!seed||!o||!output||!stats||!o->reference_uv||!o->groups||
       (o->exact_metric_xyz&&o->exact_metric_xyz!=xyz)||
       !o->group_masks||!o->group_count_guard||!o->material_guard||!o->observe||
       !o->slim.linear_solve||!o->slim.directions||o->slim.direction_weight!=16||
       o->slim.fixed_vertices||o->slim.metric_uv||o->slim.scaffold_uv||o->slim.guide_uv||
       o->slim.direction_rigid_gauge||o->slim.strip_pack||!o->slim.guard_boundary||
       o->model<0||o->model>1||o->iterations<0||o->slim.line_search<=0||o->refinement_limit<0||o->refinement_limit==INT_MAX||
       (o->refinement_limit>0&&!o->refine)||o->final_refinement_limit<0||o->final_refinement_limit==INT_MAX||o->loss_window<0||
       (o->loss_window&&(!(o->loss_relative_tolerance>0)||!isfinite(o->loss_relative_tolerance)))||
       !isfinite(o->refinement_progress_alpha)||o->refinement_progress_alpha<0||o->refinement_progress_alpha>1||
       (o->refinement_progress_alpha>0&&(!o->loss_window||!o->refinement_limit)))return -1;
    Arena_Mark mark=Arena_save(arena);volatile int status=-1;memset(stats,0,sizeof *stats);
    /* A retained feasible iterate is always exported, including later refusal. */
    memcpy(output,seed,2*nv*sizeof *output);
    TRY {
        SRSystem s;PMAMetric *reference=ARENA_ALLOC(arena,o->groups*sizeof *reference);
        PMAMetric *current_metric=ARENA_ALLOC(arena,o->groups*sizeof *current_metric),*metric=ARENA_ALLOC(arena,o->groups*sizeof *metric);
        double *current=ARENA_ALLOC(arena,2*nv*sizeof *current),*destination=ARENA_ALLOC(arena,2*nv*sizeof *destination);
        double *trial=ARENA_ALLOC(arena,2*nv*sizeof *trial),*gradient=ARENA_ALLOC(arena,2*nv*sizeof *gradient);
        SREdgeBox *boxes=NULL;memcpy(current,seed,2*nv*sizeof *current);
        if(sr_build_system(arena,xyz,nv,faces,nf,&o->slim,1,&s))RAISE(pma_failed);
        boxes=ARENA_ALLOC(arena,(s.nb?s.nb:1)*sizeof *boxes);
        if(pma_measure(&s,o->reference_uv,o,reference)||pma_measure(&s,current,o,current_metric)||current_metric[0].flips||
           sr_boundary_intersections(&s,current,boxes))RAISE(pma_failed);
        size_t failed=0,pairs=0;
        if(!pma_physical(current_metric,reference,o,&failed)||o->material_guard(o->context,current,&pairs)<=0)RAISE(pma_failed);
        stats->initial_objective=pma_objective(current_metric,o->model);stats->final_objective=stats->initial_objective;
        if(o->loss_window&&pmal_begin(stats,stats->initial_objective,current_metric[0].area))RAISE(pma_failed);
        for(int iteration=0;iteration<o->iterations;iteration++) {
            if(o->target_rms_degrees>0&&current_metric[0].mass>0&&
               sqrt(current_metric[0].angular16/(32*current_metric[0].mass))*180/3.14159265358979323846<o->target_rms_degrees) {
                stats->termination=PMA_STOP_TARGET;break;
            }
            stats->iterations=iteration+1;
            size_t old_caps=stats->refinement_limit_hits,old_failures=stats->refinement_failures,old_progress=stats->refinement_progress_stops;
            if(pma_assemble(&s,current,o))RAISE(pma_failed);
            if(o->refinement_limit>0) {
                int refined=pma_refine_step(arena,&s,current,current_metric,reference,o,boxes,iteration+1,output,stats);
                if(refined<0)RAISE(pma_failed);
                if(!refined){if(o->loss_window)stats->termination=PMA_STOP_NO_SAFE_STEP;break;}
                if(o->loss_window) {
                    int stop=pmal_update(stats,o,stats->final_objective,stats->last_proxy_derivative,
                        old_caps==stats->refinement_limit_hits&&old_failures==stats->refinement_failures&&old_progress==stats->refinement_progress_stops);
                    if(stop<0)RAISE(pma_failed);if(stop>0)break;
                }
                continue;
            }
            memcpy(destination,current,2*nv*sizeof *destination);
            if(o->slim.linear_solve(arena,o->slim.linear_context,s.nv,s.ne,s.diagonal,s.edge_a,s.edge_b,s.edge_block,s.rhs,destination))RAISE(pma_failed);
            double derivative=pma_proxy_derivative(&s,current,destination,gradient);
            double maxstep=0;for(size_t k=0;k<2*nv;k++)maxstep=fmax(maxstep,fabs(destination[k]-current[k]));
            double h=1e-5/fmax(1,maxstep),fd=0;
            for(size_t k=0;k<2*nv;k++)trial[k]=current[k]+h*(destination[k]-current[k]);
            if(pma_measure(&s,trial,o,metric))RAISE(pma_failed);fd=pma_objective(metric,o->model);
            for(size_t k=0;k<2*nv;k++)trial[k]=current[k]-h*(destination[k]-current[k]);
            if(pma_measure(&s,trial,o,metric))RAISE(pma_failed);fd=(fd-pma_objective(metric,o->model))/(2*h);
            int accepted=0;double alpha=1;
            for(int line=0;line<o->slim.line_search;line++,alpha*=.5) {
                for(size_t k=0;k<2*nv;k++)trial[k]=current[k]+alpha*(destination[k]-current[k]);
                PMAEvent event={0};event.iteration=iteration+1;event.line=line;event.alpha=alpha;
                event.previous_objective=pma_objective(current_metric,o->model);event.proxy_derivative=derivative;event.fd_derivative=fd;
                accepted=pma_evaluate(&s,current,trial,o,reference,metric,boxes,&event);stats->trials++;
                if(o->observe(o->context,trial,&event,metric,o->groups))RAISE(pma_failed);
                stats->stop_reason=event.first_failure;
                if(accepted) {
                    memcpy(current,trial,2*nv*sizeof *current);memcpy(output,current,2*nv*sizeof *output);
                    memcpy(current_metric,metric,o->groups*sizeof *metric);stats->accepted_steps++;
                    stats->final_objective=event.objective;break;
                }
            }
            fprintf(stderr,"[private metric angular] model=%d step=%d accepted=%d alpha=%.9g objective=%.12g SD=%.12g RMS=%.9g derivative=%.9g FD=%.9g stop=%d\n",
                o->model,iteration+1,accepted,accepted?alpha:0,stats->final_objective,current_metric[0].sd,
                current_metric[0].mass>0?sqrt(current_metric[0].angular16/(32*current_metric[0].mass))*180/3.14159265358979323846:NAN,
                derivative,fd,stats->stop_reason);
            if(!accepted){if(o->loss_window)stats->termination=PMA_STOP_NO_SAFE_STEP;break;}
            if(o->loss_window) {
                int stop=pmal_update(stats,o,stats->final_objective,derivative,1);
                if(stop<0)RAISE(pma_failed);if(stop>0)break;
            }
        }
        status=0;
    } EXCEPT(pma_failed) { status=-1; } FINALLY { Arena_restore(arena,mark); } END_TRY;
    if(status&&o->loss_window)stats->termination=PMA_STOP_SOLVER_FAILURE;
    stats->status=status;return status;
}
#include "closeup_angular_metric_angular_test.inc"

#include "closeup_angular_refine_test.inc"

#include "closeup_angular_loss_stop_test.inc"

#include "closeup_angular_refine_budget_test.inc"

#include "closeup_angular_global_refine_test.inc"

#include "closeup_angular_progress_refine_test.inc"
