#include "sparse_local_global.h"
/* Native private experiment: minimize the same quadratic over one displacement
 * per original source vertex. No XYZ, faces or initial cut offsets are welded.
 * H = P^T A P; b = P^T(rhs - A current); destination = current + P delta.
 * The caller's nonlinear/path/material gates retain full responsibility. */
#include "closeup_lg_metric.inc"
#include "closeup_lg_curve.inc"
#include "closeup_lg_budget.inc"
#include "closeup_lg_budget_curve.inc"
#include "closeup_lg_boundary.inc"
#include "closeup_lg_rim_learn.inc"
#include "closeup_lg_ipm_warm.inc"

typedef struct SmLinear {
    size_t nv, unique, calls;
    int32_t *map;
    SmMetric *metric;
    SmRim *rim;
    SmBudget *budget;
    double seconds, max_relative_residual;
} SmLinear;

static int sm_setup(Arena_T arena, const float *xyz, size_t nv,
                     const int32_t *aliases, int coupled, SmLinear *s)
{
    if (!arena || !xyz || !nv || nv > (size_t)INT_MAX/3 || !s) return -1;
    for (size_t v=0; v<nv; v++) {
        for (size_t k=0; k<3; k++) if (!isfinite(xyz[3*v+k])) return -1;
        if (aliases) {
            int32_t a=aliases[v];
            if (a<0 || (size_t)a>=nv || aliases[a]!=a ||
                memcmp(xyz+3*v,xyz+3*(size_t)a,3*sizeof(float))) return -1;
        }
    }
    memset(s,0,sizeof *s); s->nv=nv;
    s->map=ARENA_ALLOC(arena,nv*sizeof *s->map);
    for (size_t v=0; v<nv; v++) s->map[v]=-1;
    for (size_t v=0; v<nv; v++) if (!coupled || !aliases || aliases[v]==(int32_t)v)
        s->map[v]=(int32_t)s->unique++;
    if (coupled && aliases) for (size_t v=0; v<nv; v++) s->map[v]=s->map[aliases[v]];
    return s->unique ? 0 : -1;
}

static int sm_linear(Arena_T arena, void *context, size_t nv, size_t ne,
                      const double *diagonal, const int32_t *edge_a,
                      const int32_t *edge_b, const double *blocks,
                      const double *rhs, double *result)
{
    SmLinear *s=context;
    if (!arena || !s || nv!=s->nv || !s->map || !s->unique || !diagonal || !rhs || !result ||
        nv>(size_t)INT_MAX/3 || ne>((size_t)INT_MAX-3*nv)/4 ||
        (ne && (!edge_a || !edge_b || !blocks))) return -1;
    for (size_t v=0; v<nv; v++) {
        if (s->map[v]<0 || (size_t)s->map[v]>=s->unique) return -1;
        for (size_t k=0; k<3; k++) if (!isfinite(diagonal[3*v+k])) return -1;
        for (size_t k=0; k<2; k++) if (!isfinite(rhs[2*v+k]) || !isfinite(result[2*v+k])) return -1;
    }
    for (size_t e=0; e<ne; e++) {
        if (edge_a[e]<0 || edge_b[e]<0 || (size_t)edge_a[e]>=nv ||
            (size_t)edge_b[e]>=nv || edge_a[e]==edge_b[e]) return -1;
        for (size_t k=0; k<3; k++) if (!isfinite(blocks[3*e+k])) return -1;
    }
    if(s->metric&&s->metric->curvature&&smc_begin(s->metric->curvature,result))return -1;
    Arena_Mark mark=Arena_save(arena); volatile int rc=-1;
    double started=ves_clock_sec(); s->calls++;
    TRY {
        size_t count=3*nv+4*ne,at=0;
        int *row=ARENA_ALLOC(arena,count*sizeof *row),*col=ARENA_ALLOC(arena,count*sizeof *col);
        double *value=ARENA_ALLOC(arena,count*sizeof *value);
        double *gradient=ARENA_ALLOC(arena,2*nv*sizeof *gradient);
        double *b=ARENA_CALLOC(arena,2*s->unique,sizeof *b);
        double *delta=ARENA_CALLOC(arena,2*s->unique,sizeof *delta);
        double *residual=ARENA_ALLOC(arena,2*s->unique*sizeof *residual);
        for (size_t v=0; v<nv; v++) {
            const double *d=diagonal+3*v; int x=2*s->map[v],y=x+1;
            gradient[2*v]=rhs[2*v]-(d[0]*result[2*v]+d[1]*result[2*v+1]);
            gradient[2*v+1]=rhs[2*v+1]-(d[1]*result[2*v]+d[2]*result[2*v+1]);
            row[at]=col[at]=x; value[at++]=d[0];
            row[at]=y; col[at]=x; value[at++]=d[1];
            row[at]=col[at]=y; value[at++]=d[2];
        }
        for (size_t e=0; e<ne; e++) {
            size_t a=(size_t)edge_a[e],b0=(size_t)edge_b[e]; const double *w=blocks+3*e;
            gradient[2*a]-=w[0]*result[2*b0]+w[1]*result[2*b0+1];
            gradient[2*a+1]-=w[1]*result[2*b0]+w[2]*result[2*b0+1];
            gradient[2*b0]-=w[0]*result[2*a]+w[1]*result[2*a+1];
            gradient[2*b0+1]-=w[1]*result[2*a]+w[2]*result[2*a+1];
            int aa=s->map[a],bb=s->map[b0];
            if (aa==bb) {
                int x=2*aa,y=x+1;
                row[at]=col[at]=x; value[at++]=2*w[0];
                row[at]=y; col[at]=x; value[at++]=2*w[1];
                row[at]=col[at]=y; value[at++]=2*w[2];
            } else for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
                int x=2*aa+k,y=2*bb+l;
                row[at]=x>y?x:y; col[at]=x>y?y:x; value[at++]=w[k+l];
            }
        }
        for (size_t v=0; v<nv; v++) for (size_t k=0; k<2; k++)
            b[2*(size_t)s->map[v]+k]+=gradient[2*v+k];
        int finite=1;
        for (size_t k=0; k<at; k++) finite=finite && isfinite(value[k]);
        for (size_t k=0; k<2*s->unique; k++) finite=finite && isfinite(b[k]);
        SmBound *metric_rows=NULL;size_t bound_count=0;double *multipliers=NULL;
        if(finite&&s->metric&&smm_bounds(arena,s->metric,s->map,result,&metric_rows,&bound_count))finite=0;
        if(finite&&s->budget&&smb_bounds(arena,s->budget,s->map,result,&metric_rows,&bound_count))finite=0;
        if(finite&&s->rim&&smr_bounds(arena,s->rim,s->map,result,&metric_rows,&bound_count))finite=0;
        if(finite&&s->rim&&s->rim->learning&&smrl_bounds(arena,s->rim,s->map,result,&metric_rows,&bound_count))finite=0;
        if(finite&&s->budget&&s->budget->curvature&&smbc_begin(s->budget->curvature,result,metric_rows,bound_count))finite=0;
        if(finite&&s->metric&&s->metric->curvature&&smc_bounds(arena,s->metric,s->map,result,&metric_rows,&bound_count))finite=0;
        if(finite&&s->budget&&s->budget->curvature&&smbc_bounds(arena,s->budget,s->map,result,&metric_rows,&bound_count))finite=0;
        if(finite&&s->metric&&smo_bounds(s->metric,metric_rows,bound_count))finite=0;
        if(finite&&s->metric) {
            size_t positive[5]={0},outside[5]={0};double maximum[5]={0};
            for(size_t r=0;r<bound_count;r++) {
                int type=metric_rows[r].upper;if(type<0||type>4)continue;
                double v=metric_rows[r].rhs;positive[type]+=v>0;outside[type]+=v>1e-4+1e-12;maximum[type]=fmax(maximum[type],v);
            }
            fprintf(stderr,"[source origin] {\"call\":%zu,\"rows\":%zu,\"positive_rhs\":[%zu,%zu,%zu,%zu,%zu],\"beyond_skin\":[%zu,%zu,%zu,%zu,%zu],\"max_rhs\":[%.17g,%.17g,%.17g,%.17g,%.17g]}\n",
                s->calls,bound_count,positive[0],positive[1],positive[2],positive[3],positive[4],outside[0],outside[1],outside[2],outside[3],outside[4],maximum[0],maximum[1],maximum[2],maximum[3],maximum[4]);
        }
        if(!finite)rc=-1;
        else if(bound_count&&s->metric&&s->metric->ipm) {
            rc=smiw_project(arena,2*s->unique,at,row,col,value,b,bound_count,metric_rows,delta,&multipliers,s->metric);
        }
        else if(bound_count) rc=-1; /* no alternative constrained solver */
        else rc=Sparse_solve_sym((int)(2*s->unique),(int)at,row,col,value,b,delta,SPARSE_SPD);
        if (!rc) {
            memcpy(residual,b,2*s->unique*sizeof *b);
            for (size_t k=0; k<at; k++) {
                residual[row[k]]-=value[k]*delta[col[k]];
                if (row[k]!=col[k]) residual[col[k]]-=value[k]*delta[row[k]];
            }
            if(multipliers)for(size_t r=0;r<bound_count;r++)for(size_t k=0;k<6;k++)
                residual[metric_rows[r].index[k]]+=multipliers[r]*metric_rows[r].g[k];
            double maximum=0,scale=0;
            for (size_t k=0; k<2*s->unique; k++) {
                maximum=fmax(maximum,fabs(residual[k])); scale=fmax(scale,fabs(b[k]));
                if (!isfinite(delta[k]) || !isfinite(residual[k])) rc=-1;
            }
            double model_slope=0,quadratic=0,dual_work[5]={0};
            for(size_t k=0;k<2*s->unique;k++)model_slope-=2*b[k]*delta[k];
            for(size_t k=0;k<at;k++)quadratic+=(row[k]==col[k]?1:2)*value[k]*delta[row[k]]*delta[col[k]];
            if(multipliers)for(size_t r=0;r<bound_count;r++) {
                int type=metric_rows[r].upper;if(type>=0&&type<=4)dual_work[type]+=multipliers[r]*metric_rows[r].rhs;
            }
            fprintf(stderr,"[source origin work] {\"call\":%zu,\"model_slope\":%.17g,\"quadratic\":%.17g,\"dual_rhs_work\":[%.17g,%.17g,%.17g,%.17g,%.17g]}\n",
                s->calls,model_slope,quadratic,dual_work[0],dual_work[1],dual_work[2],dual_work[3],dual_work[4]);
            double relative=maximum/fmax(1,scale);
            s->max_relative_residual=fmax(s->max_relative_residual,relative);
            if (!isfinite(relative) || relative>1e-7) rc=-1;
            for (size_t v=0; !rc && v<nv; v++) for (size_t k=0; k<2; k++)
                if (!isfinite(result[2*v+k]+delta[2*(size_t)s->map[v]+k])) rc=-1;
            if (!rc) for (size_t v=0; v<nv; v++) for (size_t k=0; k<2; k++)
                result[2*v+k]+=delta[2*(size_t)s->map[v]+k];
        }
            if(s->metric)fprintf(stderr,"[source metric] call=%zu rows=%zu active_peak=%zu adds=%zu drops=%zu solves=%zu violation=%.9g status=%d\n",
            s->calls,bound_count,s->metric->active_peak,s->metric->adds,s->metric->drops,s->metric->solves,s->metric->max_violation,(int)rc);
            if(s->rim)fprintf(stderr,"[source boundary] call=%zu edges=%zu total_pairs=%zu max_pairs=%zu total_rows=%zu constant_rows=%zu minimum_gap=%.9g scan_seconds=%.9g status=%d\n",
            s->calls,s->rim->nb,s->rim->pairs,s->rim->max_pairs,s->rim->rows,s->rim->constant_rows,s->rim->minimum_gap,s->rim->seconds,(int)rc);
            if(s->budget)fprintf(stderr,"[source budget] call=%zu rows=%zu active_faces=%zu active_groups=%zu peak_faces=%zu peak_groups=%zu seconds=%.9g status=%d\n",
            s->calls,s->budget->rows,s->budget->active_faces,s->budget->active_groups,s->budget->peak_faces,s->budget->peak_groups,s->budget->seconds,(int)rc);
    } FINALLY {
        if(rc&&s->metric&&s->metric->ipm)s->metric->ipm->valid=0;
        Arena_restore(arena,mark);
    } END_TRY;
    s->seconds+=ves_clock_sec()-started;
    return rc;
}
