#include "asm_repair.h"
#include "../common/uv_guard.h"
#include "../common/eig3.h"
#include "../common/csr.h"
#include "../common/ves_platform.h"
#include "../common/sha256.h"
#include "../common/mesh_bin.h"
#include "asm_continuity.h"
#include "asm_flatten.h"
#include "asm_store.h"
#include "asm_report.h"
#include "asm_reading_order.h"
#include "asm_audit.h"
#include "../flatten/convex_qp.h"
#include "../flatten/active_set_qp.h"
#include "../common/pipeline_constants.h"

#include "../flatten/sparse_solve.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/ves_omp.h"

#define APR_SEAM_HUBER_LIMIT 2.0

typedef struct AprContactCopy {
    int32_t ids[6];
    int count;
    double weight;
} AprContactCopy;

typedef struct AprRank {
    size_t bad_faces,bad_seams,pairs;
    size_t unavailable_seams,fixed_seams;
    double metric_failure,seam_failure,blocking_area,area,energy;
} AprRank;

typedef struct AprSolver {
    AsmField *field;
    AsmRepairOpts opts;
    AsmRepairStats *stats;
    double *gradient, *contact_gradient, *diagonal, *face_weight;
    double *direction, *trial, *r, *z, *p, *ap;
    double *lo, *hi, *within10, *within25, *area;
    double *trial10, *trial25, *seam_scale;
    double *rigid_motion; /* chart center, translation and angle */
    uint8_t *fixed, *rigid_fixed;
    AsmSeamMeasure *seam;
    AsmContactStats contacts;
    double energy;
    size_t self_contact_regressions;
    int rigid;
    int area_clearance;
    double closure_fraction;
    AprContactCopy *local_contacts;
    size_t local_count, local_capacity;
    const double *evaluated_uv;
    const double *initial_uv; /* immutable shape alternative, never the incumbent archive */
    double separation_energy;
    int evaluated_gradient, allocation_failed;
    size_t *active_faces, n_active;
    double fixed_metric_energy;
    int metric_cached;
    SparseFactor_T factor;
    Arena_T linear_arena;
    int *free_id,n_free;
    int contact_budget_hit;   /* a contact measurement stopped at its walk budget: the solve is refused, not failed */
    double *best_uv;
    AprRank best_rank;
} AprSolver;
static int apa_initialize(AsmField *field,const uint8_t *core,uint8_t *known);
#include "asm_repair_resources.inc"
static int apr_expired(const AsmRepairOpts *o)
{return (o->deadline>0 && ves_clock_sec()>=o->deadline) || (o->workspace && o->workspace->exceeded);}
static void apr_linear_clear(AprSolver *s)
{
    Sparse_factor_free(&s->factor);
    if(s->linear_arena)Arena_dispose(&s->linear_arena);
    s->free_id=NULL;s->n_free=0;
}

void AsmRepair_defaults(AsmRepairOpts *o)
{
    memset(o,0,sizeof *o); o->iterations = 12;
    o->metric_weight = 40; o->seam_weight = 16;
    o->contact_weight = 1; o->overlap_weight = 1; o->contact_length = 2;
    o->proximal = .05; o->proximal_length = 8;
    o->regional = 1; o->admission = 1; o->clearance_iterations = 8; o->source_trials=64; o->admission_rounds=4;
    o->closure_iterations=12;
    o->context_hops=ASM_SOURCE_TRIAL_CONTEXT_HOPS;
    o->patch_faces=ASM_REPAIR_PATCH_FACES;o->patch_field_faces=ASM_REPAIR_PATCH_FIELD_FACES;
    o->patch_coordinates=ASM_REPAIR_PATCH_COORDINATES;
    o->patch_bytes=ASM_REPAIR_PATCH_BYTES;
    o->patch_seconds=45;
}
static void apr_gradient_basis(const AsmMetricFace *r, double g[6])
{
    g[2] = 1/r->length; g[3] = -r->along/(r->length*r->height);
    g[4] = 0; g[5] = 1/r->height; g[0] = -g[2]; g[1] = -g[3]-g[5];
}
/* Total overlap reduction cannot pay for a new chart self-intersection.
 * Compare each actual same-chart pair against the immutable incumbent field;
 * this also prevents exchanging existing self contacts for different pairs. */
static void apr_preserve_self_contact(void *context, size_t a, size_t b, double area)
{
    AprSolver *s=context;const AsmField *f=s->field;
    if(f->face_chart[a]!=f->face_chart[b])return;
    double old=AsmContacts_pair_faces(f->contacts,f->uv,a,b,NULL,NULL);
    if(!(old>0) || area>old+64*DBL_EPSILON*fmax(1,old))s->self_contact_regressions++;
}

/* A large unresolved seam has bounded influence. Its complete original
 * witnesses and every final audit/acceptance obligation are unchanged. */
static double apr_seam_penalty(const double r[2], double *influence)
{
    double length=hypot(r[0],r[1]),limit=APR_SEAM_HUBER_LIMIT;
    if(influence)*influence=length>limit ? limit/length : 1;
    return length>limit ? limit*(length-.5*limit) : .5*length*length;
}

#include "asm_repair_local.inc"

/* The maintained Python local/global merit: original-area ARAP, physical
 * vector Huber seams, centred triangle copies and exact intersection area.
 * Original SD/stretch distributions remain independent audit quantities. */
static double apr_objective(AprSolver *s, const double *uv, int gradient, AsmContactStats *contacts)
{
    AsmField *f = s->field; double energy = s->fixed_metric_energy; s->self_contact_regressions = 0;
    if (gradient) {
        apr_linear_clear(s);
        memset(s->gradient,0,2*f->nv*sizeof(double));
        memset(s->diagonal,0,f->nv*sizeof(double));
        for(size_t v=0;v<f->nv;v++)if(s->fixed[v])s->diagonal[v]=1;
    }
    for (size_t at = 0; at < (s->active_faces?s->n_active:f->nf); at++) {
        size_t t=s->active_faces?s->active_faces[at]:at;
        double j[4], lo, hi, det; AsmField_face(f,t,uv,j,&lo,&hi,&det);
        if (!UvGuard_positive(uv,f->faces+3*t) || !(lo > 0) || !isfinite(hi)) return INFINITY;
        double rotation[4],strain[4],square=0; apr_rotation(j,rotation);
        for(int k=0;k<4;k++){strain[k]=j[k]-rotation[k];square+=strain[k]*strain[k];}
        energy += s->opts.metric_weight*f->metric[t].area*square;
        if (!gradient) continue;
        double dj[4],g[6];
        for(int k=0;k<4;k++)dj[k]=2*s->opts.metric_weight*strain[k];
        apr_gradient_basis(f->metric+t,g);
        s->face_weight[t] = 2*s->opts.metric_weight*f->metric[t].area;
        for (int k = 0; k < 3; k++) {
            size_t v = (size_t)f->faces[3*t+k];
            s->gradient[2*v] += f->metric[t].area*(dj[0]*g[2*k]+dj[1]*g[2*k+1]);
            s->gradient[2*v+1] += f->metric[t].area*(dj[2]*g[2*k]+dj[3]*g[2*k+1]);
            s->diagonal[v] += s->face_weight[t]*(g[2*k]*g[2*k]+g[2*k+1]*g[2*k+1]);
            s->diagonal[v] += s->opts.proximal/(s->opts.proximal_length*s->opts.proximal_length)*f->metric[t].area/3;
        }
    }
    for (size_t i = 0; i < f->n_witness; i++) {
        const AsmWitness *w = f->witness+i; if (!w->valid || !w->required) continue;
        for (int d = 0; d < 2; d++) {
            double r[2]; AsmField_residual(w,d,uv,r);
            double weight = s->opts.seam_weight*w->weight, influence;
            energy += weight*apr_seam_penalty(r,&influence);
            weight *= influence;
            if (gradient) s->seam_scale[2*i+d]=weight;
            if (gradient) for (int k = 0; k < 4; k++) {
                size_t v = (size_t)w->vertex[d][k]; double c = w->coefficient[d][k];
                s->gradient[2*v] += weight*c*r[0]; s->gradient[2*v+1] += weight*c*r[1]; s->diagonal[v] += weight*c*c;
            }
        }
    }
    s->evaluated_uv=uv;s->evaluated_gradient=gradient;s->separation_energy=0;s->allocation_failed=0;
    if(gradient)s->local_count=0;
    double overlap_objective=0;
    int contact_rc=s->opts.contact_objective?
        s->opts.contact_objective(f,uv,gradient?s->contact_gradient:NULL,apr_contact_copy,s,contacts,&overlap_objective,s->opts.contact_merit_context):
        AsmContacts_objective(f->contacts,uv,gradient?s->contact_gradient:NULL,apr_contact_copy,s,contacts,&overlap_objective);
    if (contact_rc == -2) s->contact_budget_hit = 1;
    if (contact_rc || s->allocation_failed) return INFINITY;
    energy += s->opts.overlap_weight*overlap_objective+s->separation_energy;
    if (gradient) for (size_t v = 0; v < f->nv; v++) for (int d = 0; d < 2; d++) {
        s->gradient[2*v+d] += s->opts.overlap_weight*s->contact_gradient[2*v+d];
        if (s->fixed[v]) s->gradient[2*v+d] = 0;
        if (!isfinite(s->gradient[2*v+d]) || !(s->diagonal[v] > 0) || !isfinite(s->diagonal[v])) return INFINITY;
    }
    return energy;
}

static void apr_matrix(AprSolver *s, const double *x, double *y)
{
    AsmField *f = s->field;
    memset(y,0,2*f->nv*sizeof(double));
    for (size_t at = 0; at < (s->active_faces?s->n_active:f->nf); at++) {
        size_t t=s->active_faces?s->active_faces[at]:at;
        if(s->fixed[f->faces[3*t]] && s->fixed[f->faces[3*t+1]] && s->fixed[f->faces[3*t+2]])continue;
        double g[6], sum[4] = {0}; apr_gradient_basis(f->metric+t,g);
        for (int k = 0; k < 3; k++) {
            size_t v = (size_t)f->faces[3*t+k]; if (s->fixed[v]) continue;
            for (int d = 0; d < 2; d++) { sum[2*d] += g[2*k]*x[2*v+d]; sum[2*d+1] += g[2*k+1]*x[2*v+d]; }
        }
        for (int k = 0; k < 3; k++) {
            size_t v = (size_t)f->faces[3*t+k];
            for (int d = 0; d < 2; d++) y[2*v+d] += s->face_weight[t]*(g[2*k]*sum[2*d]+g[2*k+1]*sum[2*d+1]);
            double mass=s->opts.proximal/(s->opts.proximal_length*s->opts.proximal_length)*f->metric[t].area/6;
            for(int q=0;q<3;q++){
                size_t w=(size_t)f->faces[3*t+q];if(s->fixed[w])continue;
                for(int d=0;d<2;d++)y[2*v+d]+=mass*(k==q?2:1)*x[2*w+d];
            }
        }
    }
    for(size_t i=0;i<s->local_count;i++){
        const AprContactCopy *c=s->local_contacts+i;double mean[2]={0};
        for(int k=0;k<c->count;k++)if(!s->fixed[c->ids[k]])for(int d=0;d<2;d++)mean[d]+=x[2*(size_t)c->ids[k]+d]/c->count;
        for(int k=0;k<c->count;k++)for(int d=0;d<2;d++){
            size_t v=(size_t)c->ids[k];
            y[2*v+d]+=2*c->weight*((s->fixed[v]?0:x[2*v+d])-mean[d]);
        }
    }
    for (size_t i = 0; i < f->n_witness; i++) {
        const AsmWitness *w = f->witness+i; if (!w->valid || !w->required) continue;
        for (int side = 0; side < 2; side++) {
            double r[2] = {0};
            for (int k = 0; k < 4; k++) {
                size_t v = (size_t)w->vertex[side][k]; if (s->fixed[v]) continue;
                for (int d = 0; d < 2; d++) r[d] += w->coefficient[side][k]*x[2*v+d];
            }
            for (int k = 0; k < 4; k++) {
                size_t v = (size_t)w->vertex[side][k]; double wgt = s->seam_scale[2*i+side]*w->coefficient[side][k];
                for (int d = 0; d < 2; d++) y[2*v+d] += wgt*r[d];
            }
        }
    }
    for (size_t v = 0; v < f->nv; v++) if (s->fixed[v]) { y[2*v] = x[2*v]; y[2*v+1] = x[2*v+1]; }
}
static double apr_dot(const double *a, const double *b, size_t n)
{ double dot = 0; for (size_t i = 0; i < n; i++) dot += a[i]*b[i]; return dot; }
#include "asm_repair_linear.inc"

static void apr_checkpoint(AprSolver *s)
{
    AsmField *f = s->field; size_t nc = f->run->n_charts;
    memset(s->within10,0,nc*sizeof(double)); memset(s->within25,0,nc*sizeof(double)); memset(s->area,0,nc*sizeof(double));
    for (size_t t = 0; t < f->nf; t++) {
        const int32_t *face=f->faces+3*t;
        if(!s->metric_cached || !s->fixed[face[0]] || !s->fixed[face[1]] || !s->fixed[face[2]]){
            double j[4], det;AsmField_face(f,t,f->uv,j,s->lo+t,s->hi+t,&det);
        }
        size_t c = (size_t)f->face_chart[t]; double area = f->metric[t].area; s->area[c] += area;
        if (AsmMetric_within10(s->lo[t],s->hi[t])) s->within10[c] += area;
        if (AsmMetric_within25(s->lo[t],s->hi[t])) s->within25[c] += area;
    }
    s->metric_cached=1;
    for (size_t i = 0; i < f->n_seams; i++) AsmField_seam(f,i,f->uv,s->seam+i);
}

/* The two metric evaluations subtract world-sized UV coordinates before
 * dividing by source edge lengths. Their arithmetic error is not bounded
 * by a few ulps of the final singular value. For a known rigid motion,
 * propagate coordinate roundoff through both source gradient columns and
 * use the Frobenius bound on the singular-value perturbation. This applies
 * only to the non-regression comparison of an already invalid metric;
 * passing faces and the final audit use the shared small comparison band. */
static double apr_rigid_metric_roundoff(const AprSolver *s, size_t t)
{
    const AsmField *f=s->field;const AsmMetricFace *m=f->metric+t;double coordinate=1;
    for(int k=0;k<3;k++)for(int d=0;d<2;d++){
        size_t i=2*(size_t)f->faces[3*t+k]+d;
        coordinate=fmax(coordinate,fmax(fabs(f->uv[i]),fabs(s->trial[i])));
    }
    double error=32*DBL_EPSILON*coordinate;
    double first=2*error/m->length;
    double second=(2*error+fabs(m->along)*first)/m->height;
    return sqrt(2.)*hypot(first,second);
}

static const int apr_strict_band=ASM_REPAIR_STRICT_BAND;
static const int apr_masked_ends=ASM_REPAIR_MASKED_ENDS;
/* A face inside the nominal band before a step may leave it only by its own
 * arithmetic noise (see ASM_REPAIR_STRICT_BAND); faces already outside are
 * governed by the non-regression rule. */
static int apr_band_step_ok(double lo0,double hi0,double lo,double hi,double noise)
{
    if(!(lo0>=.75 && hi0<=1.25))return 1;
    return lo>=.75-noise && hi<=1.25+noise;
}

static int apr_strict_band_control(void)
{
    const double n=128*DBL_EPSILON*1.25;int failures=0;
    failures+=!apr_band_step_ok(.9,1.1,.8,1.2,n);                    /* inside stays inside */
    failures+=!apr_band_step_ok(.9,1.1,.75,1.25,n);                  /* the edge itself */
    failures+=apr_band_step_ok(.9,1.1,.8,1.25+1e-8,n);               /* the audit allowance is not a step budget */
    failures+=apr_band_step_ok(.9,1.1,.75-1e-9,1.1,n);
    failures+=!apr_band_step_ok(.9,1.1,.8,1.25+.5*n,n);              /* arithmetic noise is not a violation */
    failures+=!apr_band_step_ok(.75-5e-9,1.1,.75-6e-9,1.1,n);        /* already outside: the non-regression rule decides */
    failures+=!apr_band_step_ok(.9,1.25+5e-9,.9,1.25+5e-9,n);
    if(failures)fprintf(stderr,"  strict metric band control: %d failures\n",failures);
    return failures;
}

/* Full-domain acceptance is separate from the proxy solve. Passing original
 * faces/charts/seams keep their own budgets. Existing metric failures can only
 * improve, never borrow area from a larger neighbour. */
static int apr_guard_mode(AprSolver *s, const double *target, double alpha, int continuous)
{
    AsmField *f = s->field; size_t nc = f->run->n_charts;
    memset(s->trial10,0,nc*sizeof(double)); memset(s->trial25,0,nc*sizeof(double));
    for (size_t v = 0; v < f->nv; v++) {
        if (s->fixed[v] && (s->trial[2*v] != f->uv[2*v] || s->trial[2*v+1] != f->uv[2*v+1])) { s->stats->rejected_fixed++; return 0; }
    }
    for (size_t t = 0; t < f->nf; t++) {
        const int32_t *face=f->faces+3*t;
        double lo=s->lo[t],hi=s->hi[t];
        /* The exact coordinate check above protects every fixed vertex.
         * Reuse its already measured metric; full-sheet final auditing still
         * evaluates each original face independently. Local admission must
         * not recompute millions of unchanged SVDs for every backtrack. */
        if(!s->fixed[face[0]] || !s->fixed[face[1]] || !s->fixed[face[2]]){
        double j[4],det;AsmField_face(f,t,s->trial,j,&lo,&hi,&det);
        /* A rigid substep follows a rotation, whose determinant is one.
         * Its short endpoint chord is also checked; testing the fractional
         * chord to the full rotation would describe a different path. */
        if ((continuous && !UvGuard_interval(f->uv,s->rigid ? s->trial : target,f->faces+3*t,s->rigid ? 1 : alpha)) ||
            !UvGuard_positive(s->trial,f->faces+3*t) || !(lo > 0) || !isfinite(hi)) {
            s->stats->rejected_orientation++; return 0;
        }
        /* Share the audit's comparison band; separately bound the arithmetic
         * non-regression error of an already failing face under rigid motion. */
        double eps = fmax(ASM_METRIC_ROUNDOFF,128*DBL_EPSILON*fmax(1,fmax(hi,s->hi[t])));
        if(s->rigid && !AsmMetric_within25(s->lo[t],s->hi[t]))
            eps=fmax(eps,apr_rigid_metric_roundoff(s,t));
        double noise = 128*DBL_EPSILON*fmax(1,fmax(hi,s->hi[t]));
        if(apr_strict_band && s->rigid)noise=fmax(noise,apr_rigid_metric_roundoff(s,t));
        if ((AsmMetric_within25(s->lo[t],s->hi[t]) && !AsmMetric_within25(lo,hi)) ||
            lo+eps < fmin(.75,s->lo[t]) || hi > fmax(1.25,s->hi[t])+eps ||
            (apr_strict_band && !apr_band_step_ok(s->lo[t],s->hi[t],lo,hi,noise))) {
            s->stats->rejected_metric++;
            if (s->stats->rejected_metric == 1) {
                s->stats->first_metric_chart = f->face_chart[t];
                s->stats->first_metric_face = t-f->face_offset[f->face_chart[t]];
            }
            return 0;
        }
        }
        size_t c = (size_t)f->face_chart[t]; double area = f->metric[t].area;
        if (AsmMetric_within10(lo,hi)) s->trial10[c] += area;
        if (AsmMetric_within25(lo,hi)) s->trial25[c] += area;
    }
    for (size_t c = 0; c < nc; c++) {
        double tol = 32*DBL_EPSILON*s->area[c];
        if (s->trial10[c]+tol < fmin(.95*s->area[c],s->within10[c]) || s->trial25[c]+tol < fmin(.99*s->area[c],s->within25[c])) {
            s->stats->rejected_chart++; return 0;
        }
    }
    for (size_t i = 0; i < f->n_seams; i++) {
        if (!f->seams[i].required) continue;
        AsmSeamMeasure m; AsmField_seam(f,i,s->trial,&m); const AsmSeamMeasure *old = s->seam+i;
        /* Match GroupSeamGuard: protect a recovered complete seam. A failing
         * seam stays in the objective and audit, but its bad RMS/support is
         * not an additional immutable constraint on coupled recovery. */
        if (m.invalid != old->invalid || m.observations != old->observations || m.mass != old->mass ||
            (old->pass && !m.pass)) { s->stats->rejected_seam++; return 0; }
    }
    if(s->opts.external_guard && !s->opts.external_guard(f,s->trial,s->opts.external_context)){
        s->stats->rejected_contacts++;return 0;
    }
    return 1;
}
static int apr_guard(AprSolver *s, const double *target, double alpha)
{return apr_guard_mode(s,target,alpha,1);}

static void apr_trial(AprSolver *s, double alpha)
{
    const AsmField *f = s->field;
    if (!s->rigid) {
        for (size_t i = 0; i < 2*f->nv; i++) s->trial[i] = f->uv[i]+alpha*s->direction[i];
        return;
    }
    memcpy(s->trial,f->uv,2*f->nv*sizeof(double));
    for (size_t c = 0; c < f->run->n_charts; c++) {
        size_t start = f->vertex_offset[c];
        if (start == SIZE_MAX || s->rigid_fixed[c]) continue;
        const double *m = s->rigid_motion+5*c;
        double ct = cos(alpha*m[4]), sn = sin(alpha*m[4]);
        for (size_t k = 0; k < f->run->charts[c].nv; k++) {
            size_t v = start+k;
            double x = f->uv[2*v]-m[0], y = f->uv[2*v+1]-m[1];
            s->trial[2*v] = f->uv[2*v]+(ct-1)*x-sn*y+alpha*m[2];
            s->trial[2*v+1] = f->uv[2*v+1]+sn*x+(ct-1)*y+alpha*m[3];
        }
    }
}

static int apr_blocking_area(AprSolver *s,const double *uv,double *area)
{
    double started=ves_clock_sec();
    int rc=s->opts.contact_merit(s->field,uv,area,s->opts.contact_merit_context);
    s->stats->contact_merit_evaluations++;
    s->stats->contact_merit_seconds+=ves_clock_sec()-started;
    return rc || !isfinite(*area) || *area<0;
}

static int apr_search(AprSolver *s)
{
    AsmField *f = s->field; size_t n = 2*f->nv;
    double blocking_before=0;
    int physical_clearance=s->area_clearance && s->opts.contact_merit;
    if(physical_clearance && apr_blocking_area(s,f->uv,&blocking_before))return 0;
    /* p stores the full target and remains immutable throughout backtracking. */
    for (size_t i = 0; i < n; i++) s->p[i] = f->uv[i]+s->direction[i];
    for (int ls = 0; ls < 12 && !apr_expired(&s->opts); ls++) {
        double alpha = ldexp(1.,-ls); s->stats->trials++;
        apr_trial(s,alpha);
        if (!apr_guard(s,s->p,alpha)) continue;
        AsmContactStats contacts; double energy = apr_objective(s,s->trial,0,&contacts);
        if (!isfinite(energy)) { s->stats->rejected_energy++; continue; }
        /* Match Python's full-merit local/global acceptance. A temporary
         * interchart overlap increase can recover a source seam or metric;
         * the independent archive rank and final audit still see every pair.
         * A chart is never allowed to acquire a new self intersection. */
        if (s->self_contact_regressions) {
            s->stats->rejected_contacts++; continue;
        }
        int improved=s->area_clearance ?
            (contacts.area<s->contacts.area-1e-10 ||
             (contacts.pairs<s->contacts.pairs && contacts.area<=s->contacts.area+1e-10)) :
            energy < s->energy-1e-12*fmax(1,s->energy);
        /* Admission clearance must reduce its remaining material obligation.
         * Clearing old overlap cannot pay for worsening new material, while
         * returning within an inherited allowance need not block recovery. */
        if(physical_clearance){
            double blocking_after=0;
            if(apr_blocking_area(s,s->trial,&blocking_after)){s->stats->rejected_contacts++;continue;}
            improved=blocking_after<blocking_before-1e-12*fmax(1,blocking_before);
        }
        if (!improved) { s->stats->rejected_energy++; continue; }
        memcpy(f->uv,s->trial,n*sizeof(double)); s->energy = energy; s->contacts = contacts; return 1;
    }
    return 0;
}

static size_t apr_root(size_t *parent, size_t i)
{ while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; }
static void apr_rigid(AprSolver *s)
{
    const AsmField *f = s->field; memset(s->direction,0,2*f->nv*sizeof(double));
    memset(s->rigid_motion,0,5*f->run->n_charts*sizeof(double)); s->rigid = 1;
    for (size_t c = 0; c < f->run->n_charts; c++) {
        size_t start = f->vertex_offset[c]; if (start == SIZE_MAX || s->rigid_fixed[c]) continue;
        const AsmChart *ch = f->run->charts+c; double center[2] = {0}, mass = 0;
        for (size_t k = 0; k < ch->nv; k++) {
            size_t v = start+k; mass += f->mass[v];
            for (int d = 0; d < 2; d++) center[d] += f->mass[v]*f->uv[2*v+d];
        }
        if (!(mass > 0)) continue; center[0] /= mass; center[1] /= mass;
        double sum[2] = {0}, torque = 0, radius = 0;
        for (size_t k = 0; k < ch->nv; k++) {
            size_t v = start+k; double x = f->uv[2*v]-center[0], y = f->uv[2*v+1]-center[1];
            sum[0] += s->gradient[2*v]; sum[1] += s->gradient[2*v+1];
            torque += -y*s->gradient[2*v]+x*s->gradient[2*v+1]; radius += f->mass[v]*(x*x+y*y);
        }
        double tx = -sum[0]/mass, ty = -sum[1]/mass, angle = radius > 0 ? -torque/radius : 0;
        double length = hypot(tx,ty); if (length > 4) { tx *= 4/length; ty *= 4/length; }
        angle = fmax(-.05,fmin(.05,angle)); double ct = cos(angle), sn = sin(angle);
        double *motion = s->rigid_motion+5*c;
        motion[0] = center[0]; motion[1] = center[1]; motion[2] = tx; motion[3] = ty; motion[4] = angle;
        for (size_t k = 0; k < ch->nv; k++) {
            size_t v = start+k; double x = f->uv[2*v]-center[0], y = f->uv[2*v+1]-center[1];
            s->direction[2*v] = (ct-1)*x-sn*y+tx; s->direction[2*v+1] = sn*x+(ct-1)*y+ty;
        }
    }
}

typedef struct AprContact { AprSolver *solver; double area; size_t a,b; } AprContact;
static void apr_largest_contact(void *context, size_t a, size_t b, double area)
{
    AprContact *p = context; AsmField *f = p->solver->field;
    int32_t ca = f->face_chart[a], cb = f->face_chart[b];
    if (ca == cb || (p->solver->rigid_fixed[ca] && p->solver->rigid_fixed[cb])) return;
    if (area > p->area) { p->area = area; p->a = a; p->b = b; }
}
/* Boundary-area gradients vanish for a completely contained rigid chart.
 * An independently verified separating-axis proposal can escape that case.
 * It moves a whole chart; the same full metric/seam/contact guard still decides. */
static int apr_contact_direction(AprSolver *s)
{
    s->rigid = 0;
    AsmField *f = s->field; AprContact pair = {s,0,0,0}; AsmContactStats stats;
    if (AsmContacts_measure(f->contacts,f->uv,NULL,apr_largest_contact,&pair,&stats) || !(pair.area > 0)) return 0;
    size_t moving = pair.b, fixed = pair.a;
    if (s->rigid_fixed[f->face_chart[moving]]) { moving = pair.a; fixed = pair.b; }
    double best = DBL_MAX, dx = 0, dy = 0;
    for (int side = 0; side < 2; side++) {
        const int32_t *face = f->faces+3*(side ? moving : fixed);
        for (int k = 0; k < 3; k++) {
            size_t a = (size_t)face[k], b = (size_t)face[(k+1)%3];
            double nx = f->uv[2*b+1]-f->uv[2*a+1], ny = f->uv[2*a]-f->uv[2*b], norm = hypot(nx,ny);
            if (!(norm > 0)) continue; nx /= norm; ny /= norm;
            double alo=DBL_MAX,ahi=-DBL_MAX,blo=DBL_MAX,bhi=-DBL_MAX;
            for (int v = 0; v < 3; v++) {
                size_t va = (size_t)f->faces[3*fixed+v], vb = (size_t)f->faces[3*moving+v];
                double pa = nx*f->uv[2*va]+ny*f->uv[2*va+1], pb = nx*f->uv[2*vb]+ny*f->uv[2*vb+1];
                alo=fmin(alo,pa);ahi=fmax(ahi,pa);blo=fmin(blo,pb);bhi=fmax(bhi,pb);
            }
            double step = ahi-blo < bhi-alo ? ahi-blo+1e-6 : alo-bhi-1e-6;
            if (fabs(step) < best) { best=fabs(step);dx=step*nx;dy=step*ny; }
        }
    }
    if (!isfinite(best) || !(best > 0)) return 0;
    memset(s->direction,0,2*f->nv*sizeof(double)); size_t chart = (size_t)f->face_chart[moving], start = f->vertex_offset[chart];
    for (size_t k = 0; k < f->run->charts[chart].nv; k++) { s->direction[2*(start+k)] = dx; s->direction[2*(start+k)+1] = dy; }
    return 1;
}

#include "asm_repair_constraints.inc"

#include "asm_repair_clearance.inc"

static AprRank apr_rank(AprSolver *s)
{
    AprRank rank={0};AsmField *f=s->field;
    apr_checkpoint(s);
    for(size_t t=0;t<f->nf;t++)rank.bad_faces+=!AsmMetric_within25(s->lo[t],s->hi[t]);
    for(size_t c=0;c<f->run->n_charts;c++)if(s->area[c]>0)
        rank.metric_failure+=fmax(0,.95-s->within10[c]/s->area[c])+fmax(0,.99-s->within25[c]/s->area[c]);
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required){
        const AsmFieldSeam *seam=f->seams+i;
        /* Keep the full obligation ledger, but rank only work this field
         * can move and measure. Missing endpoints otherwise make the sum
         * infinite and hide progress on every represented incident seam. */
        if(f->vertex_offset[seam->a]==SIZE_MAX || f->vertex_offset[seam->b]==SIZE_MAX){rank.unavailable_seams++;continue;}
        if(s->opts.fixed_charts && s->opts.fixed_charts[seam->a] && s->opts.fixed_charts[seam->b]){rank.fixed_seams++;continue;}
        if(s->seam[i].pass)continue;
        const AsmSeamMeasure *m=s->seam+i;rank.bad_seams++;
        double needed=f->seams[i].original_cut?1:.9;
        if(m->invalid || !(m->mass>0) || !isfinite(m->rms))rank.seam_failure=INFINITY;
        else rank.seam_failure+=fmax(0,m->rms/ASM_SEAM_TOLERANCE_VOX-1)+fmax(0,needed-m->support/m->mass);
    }
    rank.pairs=s->contacts.pairs;rank.area=s->contacts.area;rank.energy=s->energy;
    rank.blocking_area=rank.area;
    if(s->opts.contact_merit){
        if(apr_blocking_area(s,f->uv,&rank.blocking_area))rank.blocking_area=INFINITY;
        fprintf(stderr,"[assemble checkpoint] faces %zu, metric %.9g, seams %zu, seam deficit %.9g, blocking area %.9g, total area %.9g, unavailable seams %zu, fixed seams %zu\n",
            rank.bad_faces,rank.metric_failure,rank.bad_seams,rank.seam_failure,rank.blocking_area,rank.area,rank.unavailable_seams,rank.fixed_seams);
    }
    return rank;
}
static int apr_better(AprRank a,AprRank b)
{
    if(a.bad_faces!=b.bad_faces)return a.bad_faces<b.bad_faces;
    if(a.metric_failure!=b.metric_failure)return a.metric_failure<b.metric_failure;
    if(a.bad_seams!=b.bad_seams)return a.bad_seams<b.bad_seams;
    if(a.seam_failure!=b.seam_failure)return a.seam_failure<b.seam_failure;
    if(a.blocking_area!=b.blocking_area)return a.blocking_area<b.blocking_area;
    if(a.area!=b.area)return a.area<b.area;
    return a.energy<b.energy;
}

/* Opt-in stall rule (AsmRepairOpts.stall_steps): a fitting or clearance phase
 * that keeps accepting steps without material progress on its best endpoint
 * ends there, as it would at its deadline, and the best endpoint is restored
 * as usual. Material progress: fewer bad faces or failing seams, or a
 * continuous defect (metric failure, seam failure, blocking or contact area)
 * down by stall_ratio of its value when the window opened. Contact pair
 * counts jitter with the tessellation and do not count. On PHerc0139 21^3,
 * 2,740 of the 3,648 parallel corrections that ran to their deadline had
 * made no such progress (2026-09-26). */
typedef struct AprStall { AprRank from; int steps; } AprStall;

static int apr_material(const AprRank *from,const AprRank *to,double ratio)
{
    if(to->bad_faces<from->bad_faces || to->bad_seams<from->bad_seams)return 1;
    const double a[4]={from->metric_failure,from->seam_failure,from->blocking_area,from->area};
    const double b[4]={to->metric_failure,to->seam_failure,to->blocking_area,to->area};
    for(int k=0;k<4;k++)if(a[k]>0 && b[k]<=(1-ratio)*a[k])return 1;
    return 0;
}

/* One accepted step against the phase's window; 1 when the phase stalled. */
static int apr_stall_step(AprStall *w,const AprRank *best,int steps,double ratio)
{
    if(steps<=0)return 0;
    if(apr_material(&w->from,best,ratio)){w->from=*best;w->steps=0;return 0;}
    return ++w->steps>=steps;
}

static int apr_stall_control(void)
{
    int failures=0;AprRank r={0};r.area=100;r.bad_seams=2;
    AprStall w={r,0};int stops=0,at=0;
    for(int k=1;k<=5 && !stops;k++){r.area*=.999;if(apr_stall_step(&w,&r,3,.02)){stops=1;at=k;}}
    if(!stops || at!=3)failures++;                         /* 0.1% per step: stalls on the third */
    r.area=100;w=(AprStall){r,0};stops=0;
    for(int k=1;k<=8 && !stops;k++){r.area*=.9;stops=apr_stall_step(&w,&r,3,.02);}
    if(stops)failures++;                                   /* 10% per step never stalls */
    r.area=100;r.bad_seams=2;w=(AprStall){r,0};
    apr_stall_step(&w,&r,3,.02);apr_stall_step(&w,&r,3,.02);r.bad_seams=1;
    if(apr_stall_step(&w,&r,3,.02) || w.steps)failures++;   /* a seam passing resets the window */
    w=(AprStall){r,0};if(apr_stall_step(&w,&r,0,.02))failures++;   /* off by default */
    fprintf(stderr,"  repair stall rule: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}


/* Triangle-pair counts depend on tessellation. Preserve the best physical
 * endpoint across every solve phase; a clearance phase must not erase a
 * better missing-material candidate by clearing unrelated inherited contact. */
static void apr_keep_best(AprSolver *s)
{
    if(!s->best_uv)return;
    AprRank rank=apr_rank(s);
    if(apr_better(rank,s->best_rank)){
        s->best_rank=rank;memcpy(s->best_uv,s->field->uv,2*s->field->nv*sizeof(double));
    }
}

static int apr_restore_best(AprSolver *s)
{
    if(!s->best_uv)return 0;
    memcpy(s->field->uv,s->best_uv,2*s->field->nv*sizeof(double));
    s->energy=apr_objective(s,s->field->uv,1,&s->contacts);
    return isfinite(s->energy)?0:-1;
}

#include "asm_repair_closure.inc"
#include "asm_repair_restoration.inc"
#include "asm_repair_restoration_controls.inc"

static int apr_clearance_phase(AprSolver *s,const char *phase)
{
    double phase_started=ves_clock_sec();
    const AsmRepairOpts *opts=&s->opts;s->area_clearance=1;AprStall stall={s->best_rank,0};
    if(opts->clearance_iterations && s->contacts.pairs)
        fprintf(stderr,"[assemble clearance] %s: %.9g overlap area (%zu pairs)\n",phase,s->contacts.area,s->contacts.pairs);
    for(int iteration=0;iteration<opts->clearance_iterations && s->contacts.pairs && !apr_expired(opts);iteration++){
        if(!apr_workspace_check(opts,0))break;
        apr_checkpoint(s);
        /* Contact torque can separate a long chart when opposing boundary
         * forces cancel under translation. Fit a proper rigid displacement
         * against the same source constraints before local deformation;
         * the nonlinear endpoint guard still measures every accepted move.
         * Explicit diagnostic methods continue to isolate their solver. */
        int accepted=opts->clearance_method==0 && apr_constrained(s,APC_CONTACT_POSE) && apr_search(s);
        if(!accepted && opts->clearance_method==0)accepted=apr_constrained(s,APC_AREA) && apr_search(s);
        if(!accepted && (opts->clearance_method==0 || opts->clearance_method==2))accepted=
            apr_supporting_clearance(s,fmin(1e8,1e4*pow(4,iteration)));
        if(!accepted && (opts->clearance_method==0 || opts->clearance_method==1))accepted=apr_local_clearance(s,8);
        if(!accepted && opts->clearance_method==3)accepted=apr_constrained(s,APC_AREA) && apr_search(s);
        if(!accepted)break;
        s->stats->clearance_steps++;
        s->energy=apr_objective(s,s->field->uv,1,&s->contacts);
        if(!isfinite(s->energy)){s->area_clearance=0;s->stats->clearance_seconds+=ves_clock_sec()-phase_started;return -1;}
        if(opts->finish_clearance)fprintf(stderr,"[assemble clearance] completed step %zu: %.9g area, %zu pairs, %.3f seconds\n",
            s->stats->clearance_steps,s->contacts.area,s->contacts.pairs,ves_clock_sec()-phase_started);
        apr_keep_best(s);
        if(apr_stall_step(&stall,&s->best_rank,opts->stall_steps,opts->stall_ratio)){s->stats->stall_stops++;break;}
    }
    s->area_clearance=0;s->stats->clearance_seconds+=ves_clock_sec()-phase_started;return 0;
}

static size_t apr_unresolved_incident_seams(const AprSolver *s)
{
    const AsmField *f=s->field;size_t unresolved=0;
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *seam=f->seams+i;
        if(!seam->required || !seam->complete || (s->opts.fixed_charts &&
           s->opts.fixed_charts[seam->a] && s->opts.fixed_charts[seam->b]))continue;
        /* Closure can retain an earlier iterate, so measure the retained
         * coordinates instead of relying on its last constraint checkpoint. */
        AsmSeamMeasure measure;AsmField_seam(f,i,f->uv,&measure);unresolved+=!measure.pass;
    }
    return unresolved;
}

static int apr_clearance_before_fit(AprSolver *s)
{
    if(!s->opts.clearance_first || !s->opts.clearance_iterations || !s->contacts.pairs ||
       s->opts.restoration_only || apr_expired(&s->opts))return 0;
    double started=ves_clock_sec(),allowance=ASM_REPAIR_PRE_FIT_CLEARANCE_SEC;
    double deadline=s->opts.deadline;int iterations=s->opts.clearance_iterations;
    if(deadline>0)allowance=fmin(allowance,fmax(0,deadline-started)*ASM_REPAIR_PRE_FIT_CLEARANCE_SHARE);
    s->opts.deadline=started+allowance;
    if(s->opts.clearance_iterations>ASM_REPAIR_PRE_FIT_CLEARANCE_STEPS)
        s->opts.clearance_iterations=ASM_REPAIR_PRE_FIT_CLEARANCE_STEPS;
    ConvexQp_deadline(s->opts.deadline);ActiveSetQp_deadline(s->opts.deadline);
    size_t before=s->stats->clearance_steps;int rc=apr_clearance_phase(s,"before fitting");
    s->stats->pre_fit_clearance_steps+=s->stats->clearance_steps-before;
    s->stats->pre_fit_clearance_seconds+=ves_clock_sec()-started;
    s->opts.deadline=deadline;s->opts.clearance_iterations=iterations;
    ConvexQp_deadline(deadline);ActiveSetQp_deadline(deadline);
    return rc;
}

static int apr_solve(AsmField *f, const AsmRepairOpts *opts, AsmRepairStats *st)
{
    if(!apr_workspace_check(opts,apr_workspace_estimate(f->nv,f->nf,f->n_witness,f->run->n_charts)))return 0;
    Arena_T arena = Arena_new(); AprSolver s = {0}; s.field = f; s.opts = *opts; s.stats = st;
    size_t n = 2*f->nv, nc = f->run->n_charts; int rc = -1;
    if (!f->nv || !f->nf) { Arena_dispose(&arena); return 0; }
    if(apr_expired(opts)){st->deadline_reached=1;Arena_dispose(&arena);return 0;}
    /* Phase clock for the one-line solve profile printed on exit: setup, pre-fit
     * clearance, general loop, clearance, restoration, closure, finish. */
    double mark[8];for(int k=0;k<8;k++)mark[k]=-1;mark[0]=ves_clock_sec();
    size_t walks0=0,replays0=0;AsmContacts_memo_counts(f->contacts,&walks0,&replays0);
    ConvexQp_deadline(opts->deadline);ActiveSetQp_deadline(opts->deadline);
    s.gradient=ARENA_ALLOC(arena,n*sizeof(double));s.contact_gradient=ARENA_ALLOC(arena,n*sizeof(double));
    s.direction=ARENA_ALLOC(arena,n*sizeof(double));s.trial=ARENA_ALLOC(arena,n*sizeof(double));
    s.r=ARENA_ALLOC(arena,n*sizeof(double));s.z=ARENA_ALLOC(arena,n*sizeof(double));s.p=ARENA_ALLOC(arena,n*sizeof(double));s.ap=ARENA_ALLOC(arena,n*sizeof(double));
    s.diagonal=ARENA_ALLOC(arena,f->nv*sizeof(double));s.face_weight=ARENA_ALLOC(arena,f->nf*sizeof(double));
    s.lo=ARENA_ALLOC(arena,f->nf*sizeof(double));s.hi=ARENA_ALLOC(arena,f->nf*sizeof(double));
    s.within10=ARENA_ALLOC(arena,nc*sizeof(double));s.within25=ARENA_ALLOC(arena,nc*sizeof(double));s.area=ARENA_ALLOC(arena,nc*sizeof(double));
    s.trial10=ARENA_ALLOC(arena,nc*sizeof(double));s.trial25=ARENA_ALLOC(arena,nc*sizeof(double));
    s.seam_scale=ARENA_ALLOC(arena,2*(f->n_witness ? f->n_witness : 1)*sizeof(double));
    s.rigid_motion=ARENA_ALLOC(arena,5*nc*sizeof(double));
    s.fixed=ARENA_CALLOC(arena,f->nv,1);s.rigid_fixed=ARENA_CALLOC(arena,nc,1);
    s.seam=ARENA_ALLOC(arena,(f->n_seams ? f->n_seams : 1)*sizeof(AsmSeamMeasure));
    size_t *parent=ARENA_ALLOC(arena,nc*sizeof(size_t));for(size_t c=0;c<nc;c++)parent[c]=c;
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].complete && f->seams[i].required){size_t a=apr_root(parent,f->seams[i].a),b=apr_root(parent,f->seams[i].b);parent[a>b?a:b]=a<b?a:b;}
    uint8_t *anchored=ARENA_CALLOC(arena,nc,1);
    for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX && opts->fixed_charts && opts->fixed_charts[c])anchored[apr_root(parent,c)]=1;
    for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        /* A compact field lives in the incumbent sheet's absolute frame.
         * Pinning the root of a fully included component would prevent it
         * translating away from fixed obstacles. The positive proximal mass
         * already regularizes its pose and vertex systems. */
        if(apr_root(parent,c)==c && !anchored[c] && !opts->compact){
            size_t p=f->vertex_offset[c],q=p;double longest=0;
            for(size_t k=1;k<f->run->charts[c].nv;k++){
                size_t v=p+k;double length=hypot(f->uv[2*v]-f->uv[2*p],f->uv[2*v+1]-f->uv[2*p+1]);
                if(length>longest){longest=length;q=v;}
            }
            if(q==p)goto done;
            s.fixed[p]=s.fixed[q]=1;s.rigid_fixed[c]=1;
        }
        if(opts->fixed_charts && opts->fixed_charts[c]){s.rigid_fixed[c]=1;memset(s.fixed+f->vertex_offset[c],1,f->run->charts[c].nv);}
    }
    /* Admission transports local seeds into the sheet; source transactions
     * can also move whole components. Refit alone leaves a badly interleaved
     * hierarchy after those moves, making small repairs scan the whole sheet. */
    AsmContacts_active_faces(f->contacts,NULL);
    uint8_t *moving=ARENA_CALLOC(arena,f->nf,1);size_t moving_faces=0;
    for(size_t t=0;t<f->nf;t++){for(int k=0;k<3;k++)if(!s.fixed[f->faces[3*t+k]])moving[t]=1;moving_faces+=moving[t];}
    /* A repartition is only worth its walk when a large part of the sheet
     * moves; every measurement refits the boxes anyway, so the answer is the
     * same either way (ASM_REPAIR_REBUILD_MIN_MOVING, 2026-09-16). */
    if((double)moving_faces>=ASM_REPAIR_REBUILD_MIN_MOVING*(double)f->nf && AsmContacts_rebuild(f->contacts,f->uv))goto done;
    /* ASM_REPAIR_MASKED_ENDS: with fixed charts, the context never moves, so
     * the opening and closing measurements use the solve's own mask; the
     * fixed faces' validity is still checked, face by face, below. */
    int masked_ends=apr_masked_ends && opts->fixed_charts;
    if(!masked_ends){
        s.energy=apr_objective(&s,f->uv,1,&s.contacts); if(!isfinite(s.energy)) goto done;
        st->energy_before=s.energy;st->overlap_before=s.contacts.area;st->contacts_before=s.contacts.pairs;
    }
    s.active_faces=ARENA_ALLOC(arena,f->nf*sizeof(size_t));
    for(size_t t=0;t<f->nf;t++)if(moving[t])s.active_faces[s.n_active++]=t;
    else{
        double j[4],lo,hi,det,r[4],square=0;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);apr_rotation(j,r);
        if(masked_ends && (!UvGuard_positive(f->uv,f->faces+3*t) || !(lo>0) || !isfinite(hi)))goto done;
        for(int k=0;k<4;k++)square+=(j[k]-r[k])*(j[k]-r[k]);
        s.fixed_metric_energy+=opts->metric_weight*f->metric[t].area*square;
    }
    AsmContacts_active_faces(f->contacts,moving);
    s.energy=apr_objective(&s,f->uv,1,&s.contacts);if(!isfinite(s.energy))goto done;
    if(masked_ends){st->energy_before=s.energy;st->overlap_before=s.contacts.area;st->contacts_before=s.contacts.pairs;}
    double *initial=ARENA_ALLOC(arena,n*sizeof(double));memcpy(initial,f->uv,n*sizeof(double));s.initial_uv=initial;
    double *best=ARENA_ALLOC(arena,n*sizeof(double));memcpy(best,f->uv,n*sizeof(double));
    s.best_uv=best;s.best_rank=apr_rank(&s);mark[1]=ves_clock_sec();
    /* A qualified warm start may need clearance rather than another fit.
     * Keep the resulting field as the fitting fallback's starting point;
     * the complete admission certificate still decides publication. */
    if(apr_clearance_before_fit(&s))goto done;
    mark[2]=ves_clock_sec();
    if(apr_restore_best(&s))goto done;
    if(opts->pose_iterations>0 && opts->pose_seconds>0){
        double pose_started=ves_clock_sec(),pose_stop=pose_started+opts->pose_seconds;
        for(int iteration=0;iteration<opts->pose_iterations && ves_clock_sec()<pose_stop && !apr_expired(opts);iteration++){
            if(opts->finish_clearance && !s.best_rank.bad_faces && s.best_rank.metric_failure==0 &&
               !s.best_rank.bad_seams && s.best_rank.seam_failure==0)break;
            if(!apr_workspace_check(opts,0))break;
            s.energy=apr_objective(&s,f->uv,1,&s.contacts);if(!isfinite(s.energy))goto done;
            apr_checkpoint(&s);apr_rigid(&s);
            if(!apr_constrained(&s,APC_POSE) || !apr_search(&s))break;
            st->pose_phase_steps++;st->accepted_rigid++;apr_keep_best(&s);
        }
        st->pose_phase_seconds+=ves_clock_sec()-pose_started;
        if(apr_restore_best(&s))goto done;
        fprintf(stderr,"[assemble pose] %zu joint chart-pose steps before vertex deformation, %.3f seconds\n",st->pose_phase_steps,st->pose_phase_seconds);
    }
    int general_iterations=opts->clearance_first && !s.contacts.pairs?0:opts->iterations;
    double general_started=ves_clock_sec();AprStall stall={s.best_rank,0};
    for(int iteration=0;!opts->restoration_only && iteration<general_iterations && !apr_expired(opts);iteration++){
        /* Incident-source admission seeks a feasible original patch. Once
         * its archived field has no metric, seam or contact defect, further
         * energy minimization spends the remaining cavity budget without
         * adding material. Final restoration and the complete independent
         * admission certificate still run below, including exterior contact. */
        if(opts->admission_incident_sources && !s.best_rank.bad_faces &&
           s.best_rank.metric_failure==0 && !s.best_rank.bad_seams && s.best_rank.seam_failure==0 &&
           (opts->finish_clearance || (!s.best_rank.pairs && s.best_rank.area==0 && s.best_rank.blocking_area==0))){
            fprintf(stderr,"[assemble repair] incident-source fitting stopped at %s archive after %d iterations\n",
                opts->finish_clearance?"metric/seam-feasible":"feasible",iteration);
            break;
        }
        if(!apr_workspace_check(opts,0))break;
        apr_checkpoint(&s); int accepted=0;
        apr_rigid(&s);
        if(apr_constrained(&s,APC_POSE) && apr_search(&s)){
            st->accepted_rigid++;accepted=1;
            s.energy=apr_objective(&s,f->uv,1,&s.contacts);if(!isfinite(s.energy))goto done;
            apr_checkpoint(&s);
            /* Python archives the proper-pose candidate independently. The
             * following deformation must not erase a better pose result. */
            apr_keep_best(&s);
        }
        if(!apr_expired(opts) && apr_linear(&s) && apr_constrained(&s,APC_DEFORM) && apr_search(&s)){st->accepted_deformation++;accepted=1;}
        if(!accepted && s.contacts.pairs && !apr_expired(opts) && apr_contact_direction(&s) && apr_search(&s)){st->accepted_contact++;accepted=1;}
        st->iterations++;
        fprintf(stderr,"[assemble repair] iteration %d: energy %.9g, overlap %.9g (%zu pairs), accepted %d\n",iteration+1,s.energy,s.contacts.area,s.contacts.pairs,accepted);
        fprintf(stderr,"[assemble repair] rejected trials: fixed %zu, orientation %zu, face metric %zu, chart metric %zu, seams %zu, contacts %zu, energy %zu\n",
            st->rejected_fixed,st->rejected_orientation,st->rejected_metric,st->rejected_chart,st->rejected_seam,st->rejected_contacts,st->rejected_energy);
        /* Past the deadline an iteration ends where it stands; only a
         * completed unproductive iteration is a stall. */
        if(!accepted){if(!apr_expired(opts))st->stalled=1;break;}
        s.energy=apr_objective(&s,f->uv,1,&s.contacts);if(!isfinite(s.energy))goto done;
        apr_keep_best(&s);
        if(apr_stall_step(&stall,&s.best_rank,opts->stall_steps,opts->stall_ratio)){st->stall_stops++;break;}
    }
    /* Retain the best complete encoded checkpoint, not merely the last
     * local/global iterate. Area-clearance has its own measured merit. */
    st->general_seconds+=ves_clock_sec()-general_started;mark[3]=ves_clock_sec();
    if(apr_restore_best(&s))goto done;
    if(apr_clearance_phase(&s,"before source closure"))goto done;
    mark[4]=ves_clock_sec();
    if(apr_restore_best(&s))goto done;
    if(opts->restoration_only==4){if(apr_initialize_separated(&s,opts->iterations)<0)goto done;}
    else if(opts->restoration_only==5){if(apr_restart_forest(&s)<0)goto done;}
    else if(apr_restore_source(&s))goto done;
    mark[5]=ves_clock_sec();
    apr_keep_best(&s);size_t closure_before=st->closure_steps;
    if(!opts->restoration_only && apr_close_seams(&s))goto done;
    mark[6]=ves_clock_sec();
    /* The maintained workflow clears residual contact after seam restoration.
     * Closure can change the active contacts and their available motion even
     * when it decreases overlap. Rebuild the same clearance problem on the
     * retained field once every complete incident seam passes. A partial
     * closure must remain available for continuation: the area-only clearance
     * merit does not protect a still-failing seam's recovered support.
     * An unchanged archive or a zero-iteration render needs no second pass. */
    if(opts->clearance_iterations && s.contacts.pairs && st->closure_steps>closure_before &&
       memcmp(s.best_uv,f->uv,n*sizeof(double))){
        size_t unresolved=apr_unresolved_incident_seams(&s);
        if(unresolved)fprintf(stderr,"[assemble clearance] after source closure deferred: %zu unresolved complete incident seams\n",unresolved);
        else{
            size_t clearance_before=st->clearance_steps;
            if(apr_clearance_phase(&s,"after source closure"))goto done;
            st->post_closure_clearance_steps+=st->clearance_steps-clearance_before;
        }
    }
    apr_keep_best(&s);if(apr_restore_best(&s))goto done;
    if(masked_ends)s.energy=apr_objective(&s,f->uv,0,&s.contacts);   /* still masked; usually a memo replay */
    AsmContacts_active_faces(f->contacts,NULL);s.active_faces=NULL;s.fixed_metric_energy=0;
    if(!masked_ends)s.energy=apr_objective(&s,f->uv,0,&s.contacts);
    if(!isfinite(s.energy))goto done;
    st->energy_after=s.energy;st->overlap_after=s.contacts.area;st->contacts_after=s.contacts.pairs;rc=0;
    mark[7]=ves_clock_sec();
done:
    st->deadline_reached|=apr_expired(opts);
    {   /* Where the solve's time went; '-' = phase not reached. A phase after
         * the deadline shows how long the solve kept running past it. */
        static const char *const phase[7]={"setup","prefit","general","clearance","restoration","closure","finish"};
        char text[512];size_t used=0;double last=mark[0];
        for(int k=1;k<8 && used<sizeof text;k++){
            int w=mark[k]<0?snprintf(text+used,sizeof text-used," %s -",phase[k-1]):
                  snprintf(text+used,sizeof text-used," %s %.3f",phase[k-1],mark[k]-last);
            if(w<0)break;used+=(size_t)w;if(mark[k]>=0)last=mark[k];
        }
        double over=opts->deadline>0?ves_clock_sec()-opts->deadline:0;
        size_t walks1=0,replays1=0;AsmContacts_memo_counts(f->contacts,&walks1,&replays1);
        fprintf(stderr,"[assemble solve] %zu faces, %zu pairs:%s s; %s%.3f s past the deadline; masked contact walks %zu measured, %zu replayed\n",
                f->nf,s.contacts.pairs,text,over>0?"":"not ",over>0?over:0.0,walks1-walks0,replays1-replays0);
    }
    ConvexQp_deadline(0);ActiveSetQp_deadline(0);
    if(rc && s.contact_budget_hit){
        /* 1, not -1: the field is not a sound sheet here (a stack, or a tree built for other
         * coordinates), so the caller refuses this trial and keeps the incumbent; f->uv is the
         * last accepted state (apr_search commits accepted trials only). */
        rc=1;
        fprintf(stderr,"[assemble repair] refused: the contact walk exceeded its budget (%zu leaf pairs over %zu leaves, budget %zu)\n",
                s.contacts.leaf_pairs,AsmContacts_leaves(f->contacts),AsmContacts_budget(f->contacts));
    }
    apr_linear_clear(&s);free(s.local_contacts);AsmContacts_active_faces(f->contacts,NULL);
    Arena_dispose(&arena);return rc;
}

static void apr_accumulate(AsmRepairStats *all,const AsmRepairStats *one)
{
    if(!all->regions){all->energy_before=one->energy_before;all->overlap_before=one->overlap_before;all->contacts_before=one->contacts_before;}
    all->regions++;all->energy_after=one->energy_after;all->overlap_after=one->overlap_after;all->contacts_after=one->contacts_after;
    all->accepted_rigid+=one->accepted_rigid;all->accepted_deformation+=one->accepted_deformation;all->accepted_contact+=one->accepted_contact;
    all->trials+=one->trials;all->iterations+=one->iterations;all->constrained_steps+=one->constrained_steps;all->clearance_steps+=one->clearance_steps;
    all->active_coordinates=all->active_coordinates>one->active_coordinates?all->active_coordinates:one->active_coordinates;
    all->general_seconds+=one->general_seconds;all->clearance_seconds+=one->clearance_seconds;all->factor_seconds+=one->factor_seconds;
    all->pose_phase_steps+=one->pose_phase_steps;all->pose_phase_seconds+=one->pose_phase_seconds;
    all->pre_fit_clearance_steps+=one->pre_fit_clearance_steps;all->pre_fit_clearance_seconds+=one->pre_fit_clearance_seconds;
    all->contact_merit_evaluations+=one->contact_merit_evaluations;all->contact_merit_seconds+=one->contact_merit_seconds;
    all->closure_steps+=one->closure_steps;
    all->post_closure_clearance_steps+=one->post_closure_clearance_steps;
    all->restoration_steps+=one->restoration_steps;all->restoration_trials+=one->restoration_trials;
    all->restoration_skipped_contacts+=one->restoration_skipped_contacts;
    all->restoration_levels+=one->restoration_levels;
    if(one->restoration_max_coordinates>all->restoration_max_coordinates)all->restoration_max_coordinates=one->restoration_max_coordinates;
    all->stretch_levels+=one->stretch_levels;all->stretch_steps+=one->stretch_steps;
    all->stretch_pressure_increases+=one->stretch_pressure_increases;
    all->initialization_trials+=one->initialization_trials;all->initialization_steps+=one->initialization_steps;
    all->initialization_constraints+=one->initialization_constraints;all->initialization_separated+=one->initialization_separated;
    all->forest_trials+=one->forest_trials;all->forest_initialized+=one->forest_initialized;
    all->forest_entry_restarts+=one->forest_entry_restarts;
    all->rejected_fixed+=one->rejected_fixed;all->rejected_orientation+=one->rejected_orientation;all->rejected_metric+=one->rejected_metric;
    all->rejected_chart+=one->rejected_chart;all->rejected_seam+=one->rejected_seam;all->rejected_contacts+=one->rejected_contacts;all->rejected_energy+=one->rejected_energy;
    all->stalled|=one->stalled;
    if(all->first_metric_chart<0 && one->first_metric_chart>=0){all->first_metric_chart=one->first_metric_chart;all->first_metric_face=one->first_metric_face;}
}
static void apr_source_components(const AsmField *f,size_t *parent)
{
    for(size_t c=0;c<f->run->n_charts;c++)parent[c]=c;
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && f->seams[i].complete){
        size_t a=apr_root(parent,f->seams[i].a),b=apr_root(parent,f->seams[i].b);parent[a>b?a:b]=a<b?a:b;
    }
}
static int apr_schedule(AsmField *f,const AsmRepairOpts *opts,AsmRepairStats *stats)
{
    if(!opts->regional || opts->fixed_charts){int rc=apr_solve(f,opts,stats);stats->regions=1;return rc;}
    Arena_T arena=Arena_new();size_t nc=f->run->n_charts;
    size_t *parent=ARENA_ALLOC(arena,nc*sizeof(size_t));uint8_t *fixed=ARENA_ALLOC(arena,nc);
    apr_source_components(f,parent);
    int rc=0;
    for(size_t root=0;root<nc && !apr_expired(opts);root++)if(f->vertex_offset[root]!=SIZE_MAX && apr_root(parent,root)==root){
        size_t count=0;for(size_t c=0;c<nc;c++){fixed[c]=(uint8_t)(apr_root(parent,c)!=root);if(f->vertex_offset[c]!=SIZE_MAX && !fixed[c])count++;}
        fprintf(stderr,"[assemble repair] source region %zu: %zu original charts, complete fixed-sheet collision context\n",root,count);
        AsmRepairOpts local=*opts;local.fixed_charts=fixed;local.regional=0;local.admission=0;
        AsmRepairStats part={0};part.first_metric_chart=-1;
        rc=apr_solve(f,&local,&part);apr_accumulate(stats,&part);
        if(rc==1){
            stats->regions_refused_budget++;
            fprintf(stderr,"[assemble repair] source region %zu refused by the contact walk budget; the region phase ends with the field as it stands\n",root);
            rc=0;break;
        }
        if(rc)break;
        if(opts->checkpoint && opts->checkpoint(f,opts->checkpoint_context)){rc=-1;break;}
    }
    Arena_dispose(&arena);return rc;
}

#include "asm_repair_regions.inc"
typedef struct AmrPolishStats {
    size_t rounds,vertices,trials,metric_refused,seam_refused,contact_queries;
    size_t contact_refused,query_failed,zero_gradient,pairs_before,pairs_after;
    size_t candidates,pair_increasing_steps,focused_rounds,boundary_candidates,plateau_steps;
    double area_before,area_after;
    const char *stop;
} AmrPolishStats;
static size_t amr_vertex_polish_report(AsmField *f,size_t chart,double deadline,AmrPolishStats *stats);
static size_t amr_vertex_polish(AsmField *f,size_t chart,double deadline);
static int amr_trace(AsmField *f,size_t chart,const double *before,const char *dir);
static int api_recover(AsmField *seed, AsmChart *chart, Arena_T arena,
                       const AsmRepairOpts *opts, const char **reason, FILE *trace);
static int api_cut_recover(AsmRun *trial, size_t chart, const AsmRepairOpts *opts, const char **reason, FILE *trace);
static int api_cut_recover_from(AsmRun *trial, size_t chart, const AsmRepairOpts *opts, int first, int tear, const char **reason, FILE *trace, int *prepared);
static int api_cut_stage(AsmRun *run, const AsmRun *trial, const uint8_t *core, AsmRun *staged);
#include "asm_repair_admission.inc"
#include "asm_repair_patches.inc"
#include "asm_repair_metric_seed.inc"
#include "asm_region_pack.inc"
#include "asm_metric_feasibility.inc"
#include "asm_repair_intrinsic.inc"
#include "asm_repair_metric_cut.inc"
#include "asm_repair_cut_admission.inc"
#include "asm_repair_preservation.inc"
#include "asm_repair_region_refine.inc"
static int apr_run(AsmRun *run,const AsmRepairOpts *opts,const char *dir,AsmRepairStats *stats,int defer_output);

typedef struct AprProgress { AsmRun *run;const char *dir;size_t generation; } AprProgress;
static int apr_save_progress(AsmRun *run,const char *dir)
{
    char path[2048];
    if(snprintf(path,sizeof path,"%s/stage6_progress.asr",dir)>=(int)sizeof path || AsmStore_write_run(path,run))return -1;
    Arena_T arena=Arena_new();AsmChart *display=ARENA_ALLOC(arena,run->n_charts*sizeof(AsmChart));
    memcpy(display,run->charts,run->n_charts*sizeof(AsmChart));
    for(size_t c=0;c<run->n_charts;c++){display[c].placed=AsmChart_registered(display+c);display[c].component=0;}
    snprintf(path,sizeof path,"%s/stage6_progress.png",dir);
    int rc=AsmReport_layout_png(arena,path,display,run->n_charts,NULL,0,1,8192,1);
    Arena_dispose(&arena);return rc;
}
static int apr_progress(AsmField *field,void *context)
{
    AprProgress *p=context;
    if(AsmField_commit(field,p->run) || apr_save_progress(p->run,p->dir))return -1;
    p->generation++;
    fprintf(stderr,"[assemble repair] saved recoverable checkpoint %zu: %s/stage6_progress.asr\n",p->generation,p->dir);
    return 0;
}

/* The overlap pre-check: the placed material rasterized component-blind at the
 * conflict cell (the placement's own incremental grid), and the fraction of its
 * cells that carry two charts a layer apart in 3-D.  A sheet with seams to
 * close reads ~0; a stack reads a tenth or more, and no whole-field solve can
 * make it sound. */
static double apr_stacked_fraction(const AsmRun *run,const AsmConflictOpts *o,size_t *cells,size_t *stacked)
{
    Arena_T arena=Arena_new();AsmConflictGrid *g=AsmConflictGrid_new(arena,o,run);
    for(size_t c=0;c<run->n_charts;c++)if(AsmChart_registered(run->charts+c))AsmConflictGrid_insert(g,run->charts+c);
    *stacked=AsmConflictGrid_stacked(g,cells);
    AsmConflictGrid_dispose(g);Arena_dispose(&arena);
    return *cells?(double)*stacked/(double)*cells:0.0;
}

static int apr_run(AsmRun *run,const AsmRepairOpts *opts,const char *dir,AsmRepairStats *st,int defer_output)
{
    memset(st,0,sizeof *st); st->first_metric_chart = -1;
    if(!run || !opts || opts->iterations<0 || opts->clearance_method<0 || opts->clearance_method>3 || opts->restoration_only<0 || opts->restoration_only>6 ||
        opts->compact<0 || opts->compact>1 || opts->connectivity_first<0 || opts->connectivity_first>1 ||
        !isfinite(opts->budget_seconds) || opts->budget_seconds<0 ||
        !isfinite(opts->admission_metric_margin) || opts->admission_metric_margin<0 || opts->admission_metric_margin>=.25 ||
        !isfinite(opts->patch_seconds) || opts->patch_seconds<=0 || !opts->patch_faces || opts->patch_faces>INT_MAX ||
        !opts->patch_field_faces || opts->patch_field_faces>INT_MAX || !opts->patch_coordinates || opts->patch_coordinates>INT_MAX ||
        opts->regional<0 || opts->regional>1 || opts->admission<0 || opts->admission>1 || opts->clearance_iterations<0 || opts->source_trials<0 || opts->admission_rounds<1 || opts->admission_rounds>64 || opts->closure_iterations<0 || opts->closure_iterations>256 ||
        !isfinite(opts->restoration_grid_step) || opts->restoration_grid_step<0 ||
        !(opts->metric_weight>0) || !isfinite(opts->metric_weight) || !(opts->overlap_weight>0) || !isfinite(opts->overlap_weight) ||
        !(opts->contact_length>0) || !isfinite(opts->contact_length) || !(opts->proximal_length>0) || !isfinite(opts->proximal_length) ||
        !(opts->seam_weight>0) || !isfinite(opts->seam_weight) || !(opts->contact_weight>0) || !isfinite(opts->contact_weight) ||
        !(opts->proximal>0) || !isfinite(opts->proximal))return -1;
    double start=ves_clock_sec();AsmRepairOpts working=*opts;
    if(working.budget_seconds>0 && working.deadline<=0)working.deadline=start+working.budget_seconds;
    opts=&working;AsmField *f=AsmField_new(run->arena,run);if(!f)return -1;
    /* Every objective evaluation measures the field's contacts, and that walk
     * was the whole of this stage's single-core time. */
    AsmContacts_threads(f->contacts,opts->threads);
    AprProgress progress={run,dir,0};
    if(dir && !defer_output){working.checkpoint=apr_progress;working.checkpoint_context=&progress;}opts=&working;
    st->repairable=1;
    if(opts->layer_conflict){
        double t_pre=ves_clock_sec();
        st->stacked_fraction=apr_stacked_fraction(run,opts->layer_conflict,&st->grid_cells,&st->stacked_cells);
        st->precheck_sec=ves_clock_sec()-t_pre;
        st->repairable=st->stacked_fraction<=ASM_REPAIR_MAX_OVERLAP_FRACTION;
        fprintf(stderr,"[assemble repair] overlap pre-check: %zu of %zu placed cells carry two charts a layer apart (%.4f, limit %.3f) in %.1f s: %s\n",
                st->stacked_cells,st->grid_cells,st->stacked_fraction,ASM_REPAIR_MAX_OVERLAP_FRACTION,st->precheck_sec,
                st->repairable?"repairable":"NOT REPAIRABLE -- source-repair trials and regions skipped; the sheet is the placement plus its bounded admissions");
    }
    /* Admit an already qualified source seed before context repair. Solving
     * a conflicted component here holds its unrepaired anchors fixed and can
     * spend every clearance/closure round before any context is repaired.
     * The full configured solve runs in the final admission pass, with the
     * same whole-component gates and any improved anchor field. */
    AsmRepairOpts initial=*opts;
    initial.iterations=initial.clearance_iterations=initial.closure_iterations=0;
    initial.admission_rounds=1;
    initial.layer_place=NULL;initial.layer_conflict=NULL;
    if(opts->compact){double stop=ves_clock_sec()+30;if(initial.deadline<=0 || stop<initial.deadline)initial.deadline=stop;}
    int rc=apr_admit_to(run,&initial,dir,st,"stage6_admission_initial.csv",f);
    if(!rc && st->admitted_charts){f=AsmField_new(run->arena,run);if(!f)rc=-1;}
    if(!rc && opts->compact && opts->admission && opts->iterations && dir && !apr_expired(opts)){
        /* Detect missing neighborhoods before spending the whole remaining
         * budget on contacts among already represented charts. */
        double remaining=opts->deadline>0?fmax(0,opts->deadline-ves_clock_sec()):opts->patch_seconds;
        double deadline=ves_clock_sec()+ASM_REPAIR_NEIGHBORHOOD_WORK_FRACTION*remaining;
        char path[2048];snprintf(path,sizeof path,"%s/neighborhoods",dir);ArrStats neighborhoods={0};
        rc=arr_work(run,path,deadline,ASM_REGION_REFINE_MIN_BRANCH_AREA,opts,f,&neighborhoods);
        st->admitted_charts+=neighborhoods.admitted+neighborhoods.grown;
        if(!rc && neighborhoods.admitted){f=AsmField_new(run->arena,run);if(!f)rc=-1;}
    }
    ArrHeld held={0};if(!rc && opts->compact)arr_hold(run,opts,&held);
    if(!rc && st->repairable && opts->compact)rc=app_schedule(f,opts,dir,st);
    if(!rc && st->repairable && !opts->compact)rc=apr_reconnect(f,run,opts,dir,st);
    if(!rc && st->repairable && !opts->compact)rc=apr_schedule(f,opts,st);
    if(!rc)rc=AsmField_commit(f,run);
    AsmRepairOpts final_admission=*opts;
    if(opts->compact){double stop=ves_clock_sec()+60;if(final_admission.deadline<=0 || stop<final_admission.deadline)final_admission.deadline=stop;}
    size_t before_final=st->admitted_charts;
    if(!rc)rc=apr_admit(run,&final_admission,dir,st,f);
    if(!rc && st->admitted_charts>before_final){f=AsmField_new(run->arena,run);if(!f)rc=-1;}
    if(!rc)rc=arr_refresh(run,opts,&held);arr_held_clear(&held);
    if(!rc && opts->checkpoint)rc=opts->checkpoint(f,opts->checkpoint_context);
    if(f){st->vertices=f->nv;st->faces=f->nf;st->charts=f->charts;st->obligations=f->n_seams;}
    if(!rc && dir && !defer_output)rc=AsmField_write(f,dir);
    if(dir && opts->compact){char path[2048];snprintf(path,sizeof path,"%s/stage6_budget.json",dir);FILE *fp=fopen(path,"wb");if(!fp)return -1;
        int ok=fprintf(fp,"{\"compact\":true,\"patches\":%zu,\"kept\":%zu,\"deferred\":%zu,\"oversize\":%zu,\"deadline_reached\":%s,\"moving_face_cap\":%zu,\"field_face_cap\":%zu,\"coordinate_cap\":%zu,\"patch_seconds\":%.17g,\"source_candidates\":%zu,\"source_unresolved\":%zu,\"source_unattempted\":%zu,\"seconds\":%.17g}\n",st->patches,st->patches_kept,st->patches_deferred,st->patches_oversize,st->deadline_reached?"true":"false",opts->patch_faces,opts->patch_field_faces,opts->patch_coordinates,opts->patch_seconds,st->source_candidates,st->source_unresolved,st->source_unattempted,ves_clock_sec()-start)>=0;
        if(fclose(fp))ok=0;if(!ok)rc=-1;
    }
    if(dir){char path[2048];snprintf(path,sizeof path,"%s/stage6_repair.json",dir);FILE *fp=fopen(path,"wb");if(!fp)return -1;
        char eb[64],ea[64],ob[64],oa[64],cb[64],ca[64];
        if(opts->compact){strcpy(eb,"null");strcpy(ea,"null");strcpy(ob,"null");strcpy(oa,"null");strcpy(cb,"null");strcpy(ca,"null");}
        else{
            snprintf(eb,sizeof eb,"%.17g",st->energy_before);snprintf(ea,sizeof ea,"%.17g",st->energy_after);
            snprintf(ob,sizeof ob,"%.17g",st->overlap_before);snprintf(oa,sizeof oa,"%.17g",st->overlap_after);
            snprintf(cb,sizeof cb,"%zu",st->contacts_before);snprintf(ca,sizeof ca,"%zu",st->contacts_after);
        }
        int ok=fprintf(fp,"{\"schema\":\"native-original-field-repair-v3\",\"objective\":\"original_area_arap+vector_huber_seams+centered_triangle_contact_copies+exact_pair_overlap\",\"energy_domain\":\"%s\",\"clearance_acceptance\":\"%s\",\"seam_huber_radius_vox\":%.17g,\"completed\":%s,\"vertices\":%zu,\"faces\":%zu,\"charts\":%zu,\"obligations\":%zu,\"iterations\":%zu,\"accepted_rigid\":%zu,\"accepted_deformation\":%zu,\"accepted_contact\":%zu,\"trials\":%zu,\"energy_before\":%s,\"energy_after\":%s,\"overlap_before\":%s,\"overlap_after\":%s,\"contacts_before\":%s,\"contacts_after\":%s,\"stalled\":%s,\"seconds\":%.17g,"
            "\"repairable\":%s,\"stacked_fraction\":%.6g,\"stacked_cells\":%zu,\"grid_cells\":%zu,\"precheck_seconds\":%.3f,\"regions_refused_budget\":%zu,\"trials_refused_budget\":%zu,"
            "\"rejected\":{\"fixed\":%zu,\"orientation\":%zu,\"face_metric\":%zu,\"chart_metric\":%zu,\"seam\":%zu,\"contacts\":%zu,\"energy\":%zu},\"first_metric_chart\":%d,\"first_metric_face\":%zu,"
            "\"regions\":%zu,\"constrained_steps\":%zu,\"clearance_steps\":%zu,\"source_trials\":%zu,\"recovered_connections\":%zu,\"admissions\":%zu,\"admitted_charts\":%zu,\"closure_steps\":%zu,\"post_closure_clearance_steps\":%zu,\"clearance_schedule\":\"before_source_closure_and_after_complete_changed_closure\","
            "\"restoration_steps\":%zu,\"restoration_trials\":%zu,\"restoration_skipped_contacts\":%zu,\"restoration_objective\":\"original_area_SD+physical_Huber+source_scaled_boundary\",\"restoration_coordinates\":\"bilinear_8_4_original_vertices\","
            "\"restoration_levels\":%zu,\"restoration_max_coordinates\":%zu,\"stretch_levels\":%zu,\"stretch_steps\":%zu,\"stretch_pressure_increases\":%zu,"
            "\"initialization_trials\":%zu,\"initialization_steps\":%zu,\"initialization_constraints\":%zu,\"initialization_separated\":%zu,"
            "\"forest_trials\":%zu,\"forest_initialized\":%zu,\"forest_entry_restarts\":%zu,"
            "\"settings\":{\"metric_weight\":%.17g,\"seam_weight\":%.17g,\"contact_weight\":%.17g,\"overlap_weight\":%.17g,\"contact_length\":%.17g,\"proximal\":%.17g,\"proximal_length\":%.17g,\"regional\":%d,\"admission\":%d,\"iterations\":%d,\"clearance_iterations\":%d,\"source_trials\":%d,\"admission_rounds\":%d,\"closure_iterations\":%d},\"qualification\":\"See independent stage7_audit.json\"}\n",
            opts->compact?"Per-patch fields differ; global summaries are null. See stage6_patches.csv and complete stage7_audit.json":"Placed field after initial admissions, before final admission pass",opts->compact?"Each retained patch preserves affected pair count and overlap area; no new outside contacts":"Complete affected overlap area; pair count breaks area ties",APR_SEAM_HUBER_LIMIT,rc?"false":"true",st->vertices,st->faces,st->charts,st->obligations,st->iterations,st->accepted_rigid,st->accepted_deformation,st->accepted_contact,st->trials,
            eb,ea,ob,oa,cb,ca,st->stalled?"true":"false",ves_clock_sec()-start,
            st->repairable?"true":"false",st->stacked_fraction,st->stacked_cells,st->grid_cells,st->precheck_sec,st->regions_refused_budget,st->trials_refused_budget,
            st->rejected_fixed,st->rejected_orientation,st->rejected_metric,st->rejected_chart,st->rejected_seam,st->rejected_contacts,st->rejected_energy,st->first_metric_chart,st->first_metric_face,
            st->regions,st->constrained_steps,st->clearance_steps,st->source_trials,st->recovered_connections,st->admissions,st->admitted_charts,st->closure_steps,st->post_closure_clearance_steps,
            st->restoration_steps,st->restoration_trials,st->restoration_skipped_contacts,
            st->restoration_levels,st->restoration_max_coordinates,st->stretch_levels,st->stretch_steps,st->stretch_pressure_increases,
            st->initialization_trials,st->initialization_steps,st->initialization_constraints,st->initialization_separated,
            st->forest_trials,st->forest_initialized,st->forest_entry_restarts,
            opts->metric_weight,opts->seam_weight,opts->contact_weight,opts->overlap_weight,opts->contact_length,opts->proximal,opts->proximal_length,opts->regional,opts->admission,opts->iterations,opts->clearance_iterations,opts->source_trials,opts->admission_rounds,opts->closure_iterations)>=0;
        if(fclose(fp))ok=0;if(!ok)rc=-1;
    }return rc;
}

int AsmRepair_run(AsmRun *run,const AsmRepairOpts *opts,const char *dir,AsmRepairStats *stats)
{
    return apr_run(run,opts,dir,stats,0);
}

static AsmTrim apr_test_trim(const AsmChart *chart, int32_t vertex)
{
    AsmTrim trim = {0};
    for (size_t t = 0; t < chart->nf; t++) {
        const int32_t *v = chart->faces+3*t;
        if (v[0] != vertex && v[1] != vertex && v[2] != vertex) continue;
        double a = chart->xyz[3*(size_t)v[1]]-chart->xyz[3*(size_t)v[0]], b = chart->xyz[3*(size_t)v[2]]-chart->xyz[3*(size_t)v[0]];
        double c = chart->xyz[3*(size_t)v[1]+1]-chart->xyz[3*(size_t)v[0]+1], d = chart->xyz[3*(size_t)v[2]+1]-chart->xyz[3*(size_t)v[0]+1], det = a*d-b*c;
        if (det == 0) continue;
        trim.face = (int32_t)t; trim.valid = 1; trim.inverse[0]=d/det;trim.inverse[1]=-b/det;trim.inverse[2]=-c/det;trim.inverse[3]=a/det;break;
    } return trim;
}

/* Contact regression fixtures need triangles and penetration wider than the
 * fixed material band. Scale their physical inputs, never the production policy.
 * Repeated coordinate pointers are shared by some small test charts. */
static void apr_test_scale(AsmRun *run, double scale)
{
    for(size_t c=0;c<run->n_charts;c++){
        AsmChart *ch=run->charts+c;int source=1,intrinsic=1,placed=1;
        for(size_t p=0;p<c;p++){
            source&=run->charts[p].xyz!=ch->xyz;
            intrinsic&=run->charts[p].uv!=ch->uv;
            placed&=run->charts[p].placed_uv!=ch->placed_uv;
        }
        if(source && ch->xyz)for(size_t k=0;k<3*ch->nv;k++)ch->xyz[k]=(float)(ch->xyz[k]*scale);
        if(intrinsic && ch->uv)for(size_t k=0;k<2*ch->nv;k++)ch->uv[k]=(float)(ch->uv[k]*scale);
        if(placed && ch->placed_uv)for(size_t k=0;k<2*ch->nv;k++)ch->placed_uv[k]*=scale;
        ch->pose_x*=scale;ch->pose_y*=scale;ch->area3d*=scale*scale;
    }
    for(size_t i=0;i<run->n_corr;i++){
        AsmCorr *w=run->corr+i;
        for(int k=0;k<2;k++){w->gap_a[k]=(float)(w->gap_a[k]*scale);w->gap_b[k]=(float)(w->gap_b[k]*scale);}
        for(int k=0;k<4;k++){w->trim_a.inverse[k]/=scale;w->trim_b.inverse[k]/=scale;}
    }
}

static int apr_seam_control(void)
{
    enum { N=9,NV=N*N,NF=2*(N-1)*(N-1) };
    float xyz[2][3*NV], uv[2][2*NV]; int32_t faces[3*NF]; AsmChart charts[2]={{0}};
    for (int y=0;y<N;y++) for(int x=0;x<N;x++) {
        int v=y*N+x;
        for(int side=0;side<2;side++){xyz[side][3*v]=(float)(x+10*side);xyz[side][3*v+1]=(float)y;xyz[side][3*v+2]=0;uv[side][2*v]=(float)x;uv[side][2*v+1]=(float)y;}
        if(x+1<N && y+1<N){int32_t f[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),f,sizeof f);}
    }
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.charts=charts;run.n_charts=2;
    for(int side=0;side<2;side++){AsmChart *c=charts+side;c->id=side;c->parent=-1;c->nv=NV;c->nf=NF;c->xyz=xyz[side];c->uv=uv[side];c->faces=faces;c->placed=1;c->pose_x=15*side;c->placement_state=ASM_PLACE_ROOT;}
    AsmCorr corr[N]={{0}};AsmRelation relation={0};relation.a=0;relation.b=1;relation.corr_count=N;relation.n_corr=N;relation.continuity=ASM_CONT_SOURCE;
    run.rels=&relation;run.n_rels=1;run.corr=corr;run.n_corr=N;
    for(int y=0;y<N;y++){AsmCorr *c=corr+y;c->va=y*N+N-1;c->vb=y*N;c->valid=3;c->run=0;c->gap_a[0]=c->gap_b[0]=2;
        c->trim_a=apr_test_trim(charts,c->va);c->trim_b=apr_test_trim(charts+1,c->vb);}
    int fails=0;AsmField *before=AsmField_new(arena,&run);AsmSeamMeasure original;
    if(!before){Arena_dispose(&arena);return 1;}AsmField_seam(before,0,before->uv,&original);
    if(original.mass!=8 || fabs(original.rms-5)>1e-12 || original.pass || original.observations!=N)fails++;
    AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.iterations=8;uint8_t fixed[]={1,0};opts.fixed_charts=fixed;AsmRepairStats st;
    if(AsmRepair_run(&run,&opts,NULL,&st))fails++;
    AsmField *after=AsmField_new(arena,&run);AsmSeamMeasure result;
    if(!after)fails++;
    else {
        AsmField_seam(after,0,after->uv,&result);
        if(!result.pass || result.rms>.1 || result.mass!=original.mass || !st.accepted_rigid || st.contacts_after)fails++;
        if(memcmp(before->uv,after->uv,2*NV*sizeof(double)) || memcmp(before->xyz,after->xyz,6*NV*sizeof(float)) || memcmp(before->faces,after->faces,6*NF*sizeof(int32_t)))fails++;
    }
    /* A missing direction is an unresolved retained observation, never half
     * of a passing point. A contradictory source parity also fails closed. */
    corr[0].trim_b.valid=0;AsmField *missing=AsmField_new(arena,&run);
    if(!missing)fails++;else{AsmField_seam(missing,0,missing->uv,&result);if(result.pass || result.invalid!=1 || missing->seams[0].count!=N)fails++;}
    corr[0].trim_b.valid=1;relation.flags=ASM_REL_PARITY;AsmField *parity=AsmField_new(arena,&run);
    if(!parity)fails++;else{AsmField_seam(parity,0,parity->uv,&result);if(result.pass || parity->seams[0].parity)fails++;}
    Arena_dispose(&arena);fprintf(stderr,"  native symmetric physical seam recovery: %s (%d failures)\n",fails?"FAIL":"ok",fails);return fails;
}

static int apr_rigid_backtrack_case(int conditioned)
{
    /* A compressed incumbent is allowed to rotate, but its existing minimum
     * stretch may not decrease. The midpoint of the endpoint coordinates
     * shrinks by cos(angle/2); half of the rotation preserves every edge. */
    float xyz[] = {0,0,0, 1,0,0, 0,1,0}, uv[] = {0,0, .5f,0, 0,1};
    if(conditioned){xyz[7]=.001f;uv[5]=.001f;}
    int32_t faces[] = {0,1,2};
    AsmChart chart = {0}; chart.id=0;chart.parent=-1;chart.xyz=xyz;chart.uv=uv;chart.faces=faces;
    chart.nv=3;chart.nf=1;chart.placed=1;chart.placement_state=ASM_PLACE_ROOT;
    chart.pose_x=conditioned?10000:16;chart.pose_y=conditioned?5000:0;
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.charts=&chart;run.n_charts=1;
    AsmField *f=AsmField_new(arena,&run);if(!f){Arena_dispose(&arena);return 1;}
    double trial[6],direction[6],target[6],lo[1],hi[1],within10[1],within25[1],area[1],trial10[1],trial25[1];
    double motion[]={chart.pose_x,chart.pose_y,2,-.2,.05};uint8_t fixed[3]={0},rigid_fixed[1]={0};AsmRepairStats stats={0};
    AprSolver s={0};s.field=f;s.stats=&stats;s.trial=trial;s.direction=direction;s.lo=lo;s.hi=hi;
    s.within10=within10;s.within25=within25;s.area=area;s.trial10=trial10;s.trial25=trial25;
    s.fixed=fixed;s.rigid_fixed=rigid_fixed;s.rigid_motion=motion;s.rigid=1;
    apr_checkpoint(&s);apr_trial(&s,1);memcpy(target,trial,sizeof target);
    for(int i=0;i<6;i++)direction[i]=target[i]-f->uv[i];
    int fails=0;s.rigid=0;apr_trial(&s,.5);
    if(apr_guard(&s,target,.5) || stats.rejected_metric!=1)fails++;
    s.rigid=1;
    for(int k=0;k<5;k++){
        double alpha=ldexp(1.,-k);apr_trial(&s,alpha);
        if(!apr_guard(&s,target,alpha))fails++;
        for(int a=0;a<3;a++)for(int b=a+1;b<3;b++){
            double before=hypot(f->uv[2*a]-f->uv[2*b],f->uv[2*a+1]-f->uv[2*b+1]);
            double after=hypot(trial[2*a]-trial[2*b],trial[2*a+1]-trial[2*b+1]);
            if(fabs(before-after)>(conditioned?2e-11:2e-14))fails++;
        }
    }
    /* The coordinate error bound must still reject actual contraction of
     * an already compressed face, even when a caller claims rigid motion. */
    apr_trial(&s,.5);trial[2]=.99*trial[2]+.01*trial[0];trial[3]=.99*trial[3]+.01*trial[1];
    if(apr_guard(&s,target,.5))fails++;
    rigid_fixed[0]=1;apr_trial(&s,.5);if(memcmp(trial,f->uv,sizeof trial))fails++;
    Arena_dispose(&arena);
    return fails;
}

static int apr_rigid_backtrack_control(void)
{
    int fails=apr_rigid_backtrack_case(0)+apr_rigid_backtrack_case(1);
    fprintf(stderr,"  native rigid backtracking and conditioned source metric: %s (%d failures)\n",fails?"FAIL":"ok",fails);return fails;
}

static int apr_strict_edge_control(void)
{
    int fails = 0;
    for (int upper = 0; upper < 2; upper++) for (int rigid = 0; rigid < 2; rigid++) {
        /* Roundoff at a boundary passes, but a tiny real failure cannot borrow
         * the large healthy triangle's area budget. Test both motion modes. */
        float xyz[] = {0,0,0, 1,0,0, 0,1,0, 10,0,0, 110,0,0, 10,100,0};
        float uv[] = {0,0, .75f,0, 0,1, 10,0, 110,0, 10,100};
        if (upper) uv[2] = 1.25f;
        int32_t faces[] = {0,1,2, 3,4,5};
        AsmChart chart = {0}; chart.id=0;chart.parent=-1;chart.xyz=xyz;chart.uv=uv;chart.faces=faces;
        chart.nv=6;chart.nf=2;chart.placed=1;chart.placement_state=ASM_PLACE_ROOT;
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.charts=&chart;run.n_charts=1;
        AsmField *f=AsmField_new(arena,&run);if(!f){Arena_dispose(&arena);return fails+1;}
        double trial[12],lo[2],hi[2],within10[1],within25[1],area[1],trial10[1],trial25[1];
        uint8_t fixed[6]={0};AsmRepairStats stats={0};AprSolver s={0};
        s.field=f;s.stats=&stats;s.trial=trial;s.lo=lo;s.hi=hi;s.within10=within10;s.within25=within25;
        s.area=area;s.trial10=trial10;s.trial25=trial25;s.fixed=fixed;s.rigid=rigid;
        apr_checkpoint(&s);memcpy(trial,f->uv,sizeof trial);
        if (!apr_guard(&s,trial,1)) fails++;
        for(int k=0;k<(upper ? 4 : 1);k++)trial[2]=nextafter(trial[2],upper ? INFINITY : -INFINITY);
        double j[4],minor,major,det;AsmField_face(f,0,trial,j,&minor,&major,&det);
        if (!apr_guard(&s,trial,1)) fails++;
        trial[2]=upper?1.25+2*ASM_METRIC_ROUNDOFF:.75-2*ASM_METRIC_ROUNDOFF;
        AsmField_face(f,0,trial,j,&minor,&major,&det);
        int bad=!AsmMetric_within25(minor,major), admitted=apr_guard(&s,trial,1);
        if (!bad || admitted || stats.rejected_metric!=1) fails++;
        fprintf(stderr,"  native strict metric edge: upper %d rigid %d sigma %.17g %.17g, admitted %d\n",upper,rigid,minor,major,admitted);
        Arena_dispose(&arena);
    }
    return fails;
}

typedef struct AprTestSelf { const AsmField *field; size_t pairs; double area; } AprTestSelf;
static void apr_test_self(void *context, size_t a, size_t b, double area)
{
    AprTestSelf *s=context;
    if (s->field->face_chart[a]==s->field->face_chart[b]) { s->pairs++;s->area+=area; }
}

static int apr_self_contact_control(void)
{
    /* A connected annular strip remains strictly metric-valid while its
     * endpoints curl past each other. At the same time it separates from a
     * different chart, so aggregate contact area and energy both improve. */
    enum { N=37,NV=2*N,NF=2*(N-1) };
    float xyz[2][3*NV],uv[2][2*NV];int32_t faces[3*NF];AsmChart charts[2]={{0}};
    double target[4*NV];
    for (int i=0;i<N;i++) for (int side=0;side<2;side++) {
        int v=2*i+side;double angle=(300*3.141592653589793/180)*i/(N-1),radius=10+side;
        for(int c=0;c<2;c++) {
            xyz[c][3*v]=(float)(radius*cos(angle));xyz[c][3*v+1]=(float)(radius*sin(angle));xyz[c][3*v+2]=(float)c;
            uv[c][2*v]=(float)(.85*xyz[c][3*v]);uv[c][2*v+1]=(float)(.85*xyz[c][3*v+1]);
        }
        double next_angle=(370*3.141592653589793/180)*i/(N-1);
        target[2*v]=.5+.8*radius*cos(next_angle);target[2*v+1]=.8*radius*sin(next_angle);
        if(i+1<N && side==0){int32_t tri[]={v,v+1,v+3,v,v+3,v+2};memcpy(faces+6*i,tri,sizeof tri);}
    }
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.charts=charts;run.n_charts=2;
    for(int c=0;c<2;c++){AsmChart *ch=charts+c;ch->id=c;ch->parent=-1;ch->nv=NV;ch->nf=NF;ch->xyz=xyz[c];ch->uv=uv[c];ch->faces=faces;ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;}
    apr_test_scale(&run,64);for(int k=0;k<2*NV;k++)target[k]*=64;
    AsmField *f=AsmField_new(arena,&run);if(!f){Arena_dispose(&arena);return 1;}
    memcpy(target+2*NV,f->uv+2*NV,2*NV*sizeof(double));target[0]=f->uv[0];target[1]=f->uv[1];
    AprSolver s={0};AsmRepairStats stats={0};s.field=f;s.stats=&stats;AsmRepair_defaults(&s.opts);
    s.trial=ARENA_ALLOC(arena,4*NV*sizeof(double));s.direction=ARENA_ALLOC(arena,4*NV*sizeof(double));s.p=ARENA_ALLOC(arena,4*NV*sizeof(double));
    s.lo=ARENA_ALLOC(arena,2*NF*sizeof(double));s.hi=ARENA_ALLOC(arena,2*NF*sizeof(double));
    s.within10=ARENA_ALLOC(arena,2*sizeof(double));s.within25=ARENA_ALLOC(arena,2*sizeof(double));s.area=ARENA_ALLOC(arena,2*sizeof(double));
    s.trial10=ARENA_ALLOC(arena,2*sizeof(double));s.trial25=ARENA_ALLOC(arena,2*sizeof(double));s.fixed=ARENA_CALLOC(arena,2*NV,1);
    s.fixed[0]=1;memset(s.fixed+NV,1,NV);
    apr_checkpoint(&s);s.energy=apr_objective(&s,f->uv,0,&s.contacts);
    AsmContactStats before=s.contacts,proposed;AprTestSelf original={f,0,0},wanted={f,0,0};
    AsmContacts_measure(f->contacts,f->uv,NULL,apr_test_self,&original,&before);
    double next_energy=apr_objective(&s,target,0,&proposed);
    AsmContacts_measure(f->contacts,target,NULL,apr_test_self,&wanted,&proposed);
    memcpy(s.trial,target,4*NV*sizeof(double));int shape=apr_guard(&s,target,1);
    int fails=original.pairs || !wanted.pairs || !(proposed.area<before.area) || !(next_energy<s.energy) || !shape;
    for(int i=0;i<4*NV;i++)s.direction[i]=target[i]-f->uv[i];
    int accepted=apr_search(&s);AsmContactStats after;AprTestSelf final={f,0,0};
    AsmContacts_measure(f->contacts,f->uv,NULL,apr_test_self,&final,&after);
    if(final.pairs || (accepted && !(after.area<before.area)))fails++;
    fprintf(stderr,"  native self-contact exchange: metric-valid %d, overlap %.9g -> proposed %.9g, new self %zu, accepted %d final self %zu, failures %d\n",shape,before.area,proposed.area,wanted.pairs,accepted,final.pairs,fails);
    Arena_dispose(&arena);return fails;
}
static int apr_seam_penalty_selftest(void)
{
    const double controls[][2]={{0,0},{1,.5},{2,0},{3,4},{-30,40}};int fails=0;
    const double expected[]={0,.625,2,8,98};
    for(size_t k=0;k<sizeof controls/sizeof *controls;k++){
        double influence,energy=apr_seam_penalty(controls[k],&influence);
        if(fabs(energy-expected[k])>1e-12 || !(influence>0) || influence>1)fails++;
        double rotated[2]={-controls[k][1],controls[k][0]},rotated_scale;
        if(apr_seam_penalty(rotated,&rotated_scale)!=energy || rotated_scale!=influence)fails++;
        if(hypot(influence*controls[k][0],influence*controls[k][1])>APR_SEAM_HUBER_LIMIT+1e-12)fails++;
        for(int d=0;d<2;d++){
            double r[2]={controls[k][0],controls[k][1]},step=1e-6;
            r[d]+=step;double plus=apr_seam_penalty(r,NULL);
            r[d]-=2*step;double minus=apr_seam_penalty(r,NULL);
            if(fabs((plus-minus)/(2*step)-influence*controls[k][d])>2e-6)fails++;
        }
    }
    fprintf(stderr,"  native robust seam objective and gradient: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}

/* Shared fixture for native repair and solver concurrency controls. */
static int ags_fixture(Arena_T arena,AsmRun *run);

#include "asm_repair_fixture_test.inc"
#include "asm_repair_region_refine_test.inc"
#include "asm_repair_admission_test.inc"
#include "asm_repair_patch_test.inc"
#include "asm_region_pack_test.inc"
#include "asm_sheet_review.inc"
#include "asm_repair_joint_frontier_test.inc"
#include "asm_repair_rank_test.inc"
#include "asm_repair_failure_cache_test.inc"
#include "asm_repair_metric_seed_test.inc"
#include "asm_repair_metric_cut_test.inc"
#include "asm_repair_intrinsic_test.inc"
#include "asm_repair_cut_admission_test.inc"

#include "asm_repair_resources_test.inc"

int AsmRepair_selftest(void)
{
    /* A curved, isometric source grid with a smooth non-rigid UV distortion.
     * The solve must improve the original metric without deleting faces,
     * changing source bytes, folding, or moving an explicitly frozen chart. */
    enum { N=9, NV=N*N, NF=2*(N-1)*(N-1) };
    float xyz[3*NV], uv[2*NV]; int32_t faces[3*NF];
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){int v=y*N+x;double angle=.08*x;
        xyz[3*v]=(float)(sin(angle)/.08);xyz[3*v+1]=(float)y;xyz[3*v+2]=(float)((1-cos(angle))/.08);
        uv[2*v]=(float)(x+.14*sin(3.141592653589793*x/(N-1))*sin(3.141592653589793*y/(N-1)));uv[2*v+1]=(float)y;
        if(x+1<N && y+1<N){int at=6*(y*(N-1)+x);int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+at,tri,sizeof tri);}}
    Arena_T arena=Arena_new();AsmChart chart={0};chart.id=0;chart.parent=-1;chart.nv=NV;chart.nf=NF;chart.xyz=xyz;chart.uv=uv;chart.faces=faces;chart.placed=1;chart.placement_state=ASM_PLACE_ROOT;
    AsmRun run={0};run.arena=arena;run.charts=&chart;run.n_charts=1;
    AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.iterations=5;AsmRepairStats st;int fails=ConvexQp_selftest()+ActiveSetQp_selftest()+(Sparse_selftest()!=0)+apr_factor_failure_control()+apr_causal_control()+apr_pose_model_control()+apr_strict_edge_control()+apr_self_contact_control()+apr_seam_penalty_selftest()+apr_rank_control()+apr_admission_control()+apr_admission_locality_control()+apr_solve_determinism_control()+apa_batch_cap_control()+apr_source_reconstruction_control();
    fails+=asv_selftest()+apr_frontier_admission_control()+apr_overlay_contact_control()+apr_overlay_replacement_control()+apr_overlay_incremental_control()+apr_overlay_query_budget_control()+apr_overlay_threads_control()+apr_frontier_obligation_control()+anp_selftest()+ajf_selftest()+apa_start_area_control()+apr_strict_band_control()+arr_selftest()+ags_solver_owner_control()+afc_selftest()+amb_selftest()+amf_selftest()+acut_selftest()+api_selftest()+api_cut_selftest();
    float saved_xyz[3*NV];int32_t saved_faces[3*NF];memcpy(saved_xyz,xyz,sizeof xyz);memcpy(saved_faces,faces,sizeof faces);
    if(AsmRepair_run(&run,&opts,NULL,&st) || !chart.placed_uv || !(st.energy_after<.95*st.energy_before) || !st.accepted_deformation ||
        st.contacts_after || memcmp(xyz,saved_xyz,sizeof xyz) || memcmp(faces,saved_faces,sizeof faces))fails++;
    double saved_uv[2*NV];if(chart.placed_uv){memcpy(saved_uv,chart.placed_uv,sizeof saved_uv);uint8_t fixed[]={1};opts.fixed_charts=fixed;
        if(AsmRepair_run(&run,&opts,NULL,&st) || memcmp(saved_uv,chart.placed_uv,sizeof saved_uv))fails++;}
    Arena_dispose(&arena);fails+=apr_workspace_selftest();fails+=app_selftest();fails+=apr_seam_control();fails+=apr_rigid_backtrack_control();fails+=apr_admission_forest_control();fails+=apr_closure_plane_control();fails+=apr_closure_intrinsic_control();
    fails+=ars_derivative_control();fails+=ars_boundary_control();fails+=ars_stretch_control();fails+=ars_entry_forest_control();fails+=apr_stall_control();
    fprintf(stderr,"  native coupled repair: %s (%d failures)\n",fails?"FAIL":"ok",fails);return fails;
}
