/* A local possibility certificate for avoiding a second cone ray. The
 * temporary fan is discarded; this cannot certify the complete material. */
static int ic_slim_linear(Arena_T arena, void *context, size_t nv, size_t ne,
                          const double *diagonal, const int32_t *edge_a,
                          const int32_t *edge_b, const double *blocks,
                          const double *rhs, double *result);

/* Identical bank opening for the complete cut and its temporary fan checks. */
static int ic_cone_open_banks(Arena_T arena, MeshBinData *m,
                              const int32_t *mapping, size_t source_nv,
                              double *alpha)
{
    Arena_Mark mark=Arena_save(arena);
    size_t nv=m->nv,nf=m->nf;
    double *normals=ARENA_CALLOC(arena,3*nv,sizeof *normals);
    double *displacement=ARENA_CALLOC(arena,2*nv,sizeof *displacement);
    uint32_t *copies=ARENA_CALLOC(arena,source_nv,sizeof *copies);
    uint64_t *keys=ARENA_ALLOC(arena,3*nf*sizeof *keys);
    for (size_t v=0;v<nv;v++) copies[mapping[v]]++;
    for (size_t v=0;v<nv;v++) normals[3*v+2]=INFINITY;
    for (size_t f=0;f<nf;f++) for (int k=0;k<3;k++) {
        uint32_t a=(uint32_t)m->faces[3*f+k],b=(uint32_t)m->faces[3*f+(k+1)%3];
        keys[3*f+k]=((uint64_t)(a<b?a:b)<<32)|(uint64_t)(a<b?b:a);
    }
    qsort(keys,3*nf,sizeof *keys,ic_cone_key_compare);
    size_t nb=0;
    for (size_t k=0;k<3*nf;) {
        size_t end=k+1;while (end<3*nf && keys[end]==keys[k]) end++;
        if (end-k==1) keys[nb++]=keys[k];
        k=end;
    }
    for (size_t f=0;f<nf;f++) for (int k=0;k<3;k++) {
        size_t a=(size_t)m->faces[3*f+k],b=(size_t)m->faces[3*f+(k+1)%3],c=(size_t)m->faces[3*f+(k+2)%3];
        uint64_t key=((uint64_t)(a<b?a:b)<<32)|(uint64_t)(a<b?b:a);
        if (!ic_cone_key_has(keys,nb,key)) continue;
        double x=m->uv64[2*b]-m->uv64[2*a],y=m->uv64[2*b+1]-m->uv64[2*a+1],length=hypot(x,y);
        double det=x*(m->uv64[2*c+1]-m->uv64[2*a+1])-y*(m->uv64[2*c]-m->uv64[2*a]);
        double orientation=det>0?1:-1;
        if (!(length>0)) continue;
        for (int endpoint=0;endpoint<2;endpoint++) {
            size_t v=endpoint?b:a;
            normals[3*v]+=-orientation*y/length;normals[3*v+1]+=orientation*x/length;
            normals[3*v+2]=fmin(normals[3*v+2],length);
        }
    }
    for (size_t v=0;v<nv;v++) if (copies[mapping[v]]>1) {
        double length=hypot(normals[3*v],normals[3*v+1]);
        if (length>0 && isfinite(normals[3*v+2])) for (int q=0;q<2;q++)
            displacement[2*v+q]=.02*normals[3*v+2]*normals[3*v+q]/length;
    }
    int rc=SlimRefine_open_seam_seed(arena,m->verts,nv,m->faces,nf,m->uv64,displacement,m->uv64,alpha);
    Arena_restore(arena,mark);return rc;
}

enum { IC_FAN_PASS,IC_FAN_SUPPORT,IC_FAN_SEED,IC_FAN_SOLVE,IC_FAN_GUARD,
       IC_FAN_METRIC,IC_FAN_DIRECTION,IC_FAN_REASONS };
typedef struct {
    Arena_T arena;
    const MeshBinData *mesh;
    const SlimRefineDirection *directions;
    size_t reasons[IC_FAN_REASONS];
    size_t guard_reverted,guard_metric,guard_objective,guard_embedding;
    size_t seed_certificates,optimized_certificates;
    size_t weight_certificates[CLOSEUP_CONE_FAN_WEIGHT_STAGES];
    size_t metric_refused_with_direction_pass;
    double min_sigma,max_sigma,max_rms;
} IcConeFanContext;

/* UV Jacobian relative to each unchanged physical triangle, with no fitted
 * rest metric. The reciprocal determinant form avoids cancellation at small
 * minimum singular values. */
static int ic_cone_fan_metric(const MeshBinData *m,const double *uv,
                               double *minimum,double *maximum)
{
    *minimum=INFINITY;*maximum=0;
    for (size_t f=0;f<m->nf;f++) {
        size_t a=(size_t)m->faces[3*f],b=(size_t)m->faces[3*f+1],c=(size_t)m->faces[3*f+2];
        double x[3],y[3],cross[3],xx=0,xy=0;
        for (int k=0;k<3;k++) { x[k]=(double)m->verts[3*b+k]-m->verts[3*a+k];y[k]=(double)m->verts[3*c+k]-m->verts[3*a+k];xx+=x[k]*x[k];xy+=x[k]*y[k]; }
        for (int k=0;k<3;k++) cross[k]=x[(k+1)%3]*y[(k+2)%3]-x[(k+2)%3]*y[(k+1)%3];
        double length=sqrt(xx),height=sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2])/length;
        if (!(length>0 && height>0) || !isfinite(height)) return -1;
        double j[4]={(uv[2*b]-uv[2*a])/length,0,(uv[2*b+1]-uv[2*a+1])/length,0};
        j[1]=(uv[2*c]-uv[2*a]-j[0]*xy/length)/height;
        j[3]=(uv[2*c+1]-uv[2*a+1]-j[2]*xy/length)/height;
        double aa=j[0]*j[0]+j[2]*j[2],bb=j[1]*j[1]+j[3]*j[3],ab=j[0]*j[1]+j[2]*j[3];
        double largest=sqrt(.5*(aa+bb+hypot(aa-bb,2*ab)));
        double smallest=fabs(j[0]*j[3]-j[1]*j[2])/largest;
        if (!(smallest>0) || !isfinite(largest)) return -1;
        *minimum=fmin(*minimum,smallest);*maximum=fmax(*maximum,largest);
    }
    return 0;
}

static int ic_cone_fan_certificate(const MeshBinData *fan,const double *uv,
                                    const SlimRefineStats *stats,
                                    double direction_weight,
                                    double *minimum,double *maximum,double *rms)
{
    double reference_scale=CLOSEUP_FIBER_DIRECTION_WEIGHT/direction_weight;
    double reference_before=stats->energy_before+reference_scale*stats->direction_energy_before;
    double reference_after=stats->energy_after+reference_scale*stats->direction_energy_after;
    if (stats->reverted || stats->flips_after || stats->boundary_intersections_after ||
        !isfinite(stats->energy_after) || stats->energy_after>stats->energy_before+1e-8 ||
        !isfinite(stats->objective_after) || stats->objective_after>stats->objective_before+1e-8 ||
        !isfinite(reference_after) || reference_after>reference_before+1e-8) return IC_FAN_GUARD;
    if (stats->direction_faces!=fan->nf || stats->direction_area_fraction<1-1e-12) return IC_FAN_SUPPORT;
    *rms=sqrt(stats->direction_angle_after)*180/3.14159265358979323846;
    if (ic_cone_fan_metric(fan,uv,minimum,maximum) ||
        *minimum<1/CLOSEUP_CONE_FAN_SCALE || *maximum>CLOSEUP_CONE_FAN_SCALE) return IC_FAN_METRIC;
    if (!isfinite(*rms) || *rms>CLOSEUP_FIBER_RMS_DEG) return IC_FAN_DIRECTION;
    return IC_FAN_PASS;
}

static int ic_cone_fan_feasible(void *context,int32_t vertex,
                                const int32_t *corners,const int32_t *rim,size_t nf)
{
    IcConeFanContext *check=context;
    const MeshBinData *source=check->mesh;
    Arena_T arena=check->arena;
    Arena_Mark mark=Arena_save(arena);
    int reason=IC_FAN_SUPPORT;
    double minimum=0,maximum=0,rms=INFINITY;
    if (!check->directions || (!source->uv && !source->uv64)) goto done;
    for (size_t f=0;f<nf;f++) {
        const SlimRefineDirection *d=check->directions+(size_t)corners[f]/3;
        if (!(d->weight>0) || !isfinite(d->weight)) goto done;
        for (int k=0;k<3;k++) if (!isfinite(d->h[k]) || !isfinite(d->v[k])) goto done;
    }
    MeshBinData fan={0};fan.nv=nf+2;fan.nf=nf;
    fan.verts=ARENA_ALLOC(arena,3*fan.nv*sizeof *fan.verts);
    fan.faces=ARENA_ALLOC(arena,3*nf*sizeof *fan.faces);
    fan.uv64=ARENA_ALLOC(arena,2*fan.nv*sizeof *fan.uv64);
    int32_t *mapping=ARENA_ALLOC(arena,fan.nv*sizeof *mapping);
    int32_t *local_alias=ARENA_ALLOC(arena,fan.nv*sizeof *local_alias);
    SlimRefineDirection *directions=ARENA_ALLOC(arena,nf*sizeof *directions);
    double *result=ARENA_ALLOC(arena,2*fan.nv*sizeof *result);
    mapping[0]=vertex;
    for (size_t v=1;v<fan.nv;v++) mapping[v]=rim[v-1];
    for (size_t v=0;v<fan.nv;v++) {
        size_t original=(size_t)mapping[v];
        memcpy(fan.verts+3*v,source->verts+3*original,3*sizeof(float));
        for (int k=0;k<2;k++) fan.uv64[2*v+k]=MeshBin_uv(source,2*original+k)-MeshBin_uv(source,2*(size_t)vertex+k);
        /* Compact source aliases: avoid clearing source.nv entries per fan. */
        local_alias[v]=(int32_t)v;
        for (size_t p=0;p<v;p++) if (mapping[p]==mapping[v]) { local_alias[v]=local_alias[p];break; }
    }
    for (size_t f=0;f<nf;f++) {
        size_t original=(size_t)corners[f]/3;
        for (int k=0;k<3;k++) {
            int32_t v=source->faces[3*original+k];
            if (v==vertex) fan.faces[3*f+k]=0;
            else if (v==rim[f]) fan.faces[3*f+k]=(int32_t)f+1;
            else if (v==rim[f+1]) fan.faces[3*f+k]=(int32_t)f+2;
            else { reason=IC_FAN_GUARD;goto done; }
        }
        directions[f]=check->directions[original];
    }
    reason=IC_FAN_SEED;
    double alpha=0;
    if (ic_cone_open_banks(arena,&fan,local_alias,fan.nv,&alpha)) goto done;
    SlimRefineOpts opts;SlimRefineStats stats={0};SlimRefine_defaults(&opts);
    opts.energy_tolerance=CLOSEUP_FLAT_TOL;
    opts.strip_pack=0;
    opts.local_steps=1;opts.local_step_min_faces=1;
    opts.local_relaxation=4;opts.local_relax_boundary=1;
    opts.linear_solve=ic_slim_linear;
    opts.directions=directions;opts.direction_weight=CLOSEUP_FIBER_DIRECTION_WEIGHT;
    opts.direction_rigid_gauge=1;
    /* First certify the existing opened fan after an exact rigid gauge. A
     * later joint objective step can fail pure-SD nonworsening even though
     * this independently measured seed already satisfies both strict gates. */
    for (int optimize=0;optimize<=CLOSEUP_CONE_FAN_WEIGHT_STAGES;optimize++) {
        opts.iterations=optimize?CLOSEUP_CONE_FAN_ITERS:0;
        opts.direction_weight=ldexp(CLOSEUP_FIBER_DIRECTION_WEIGHT,optimize?1-optimize:0);
        reason=IC_FAN_SOLVE;
        if (SlimRefine_run_double(arena,fan.verts,fan.nv,fan.faces,nf,fan.uv64,&opts,result,&stats)) break;
        reason=ic_cone_fan_certificate(&fan,result,&stats,opts.direction_weight,&minimum,&maximum,&rms);
        if (reason==IC_FAN_PASS) {
            if (optimize) { check->optimized_certificates++;check->weight_certificates[optimize-1]++; }
            else check->seed_certificates++;
            break;
        }
        if (reason==IC_FAN_SUPPORT) break;
    }
    if (reason==IC_FAN_GUARD) {
        check->guard_reverted+=stats.reverted!=0;
        check->guard_metric+=!isfinite(stats.energy_after) || stats.energy_after>stats.energy_before+1e-8;
        double reference_scale=CLOSEUP_FIBER_DIRECTION_WEIGHT/opts.direction_weight;
        double reference_before=stats.energy_before+reference_scale*stats.direction_energy_before;
        double reference_after=stats.energy_after+reference_scale*stats.direction_energy_after;
        check->guard_objective+=!isfinite(stats.objective_after) || stats.objective_after>stats.objective_before+1e-8 ||
                                !isfinite(reference_after) || reference_after>reference_before+1e-8;
        check->guard_embedding+=stats.flips_after!=0 || stats.boundary_intersections_after!=0;
    }
done:
    check->reasons[reason]++;
    if (reason==IC_FAN_METRIC && isfinite(rms) && rms<=CLOSEUP_FIBER_RMS_DEG) check->metric_refused_with_direction_pass++;
    if (reason==IC_FAN_PASS) {
        if (check->reasons[IC_FAN_PASS]==1 || minimum<check->min_sigma) check->min_sigma=minimum;
        check->max_sigma=fmax(check->max_sigma,maximum);check->max_rms=fmax(check->max_rms,rms);
    }
    Arena_restore(arena,mark);return reason==IC_FAN_PASS;
}

static void ic_cone_fan_report(const IcConeFanContext *check)
{
    size_t tested=0;for (int k=0;k<IC_FAN_REASONS;k++) tested+=check->reasons[k];
    fprintf(stderr,"[ink_closeup] single-slit certificates: tested=%zu accepted=%zu "
            "refused-support=%zu seed=%zu solve=%zu guard=%zu metric=%zu direction=%zu; "
            "accepted-sigma=[%.9g,%.9g] max-rms-deg=%.9g seed=%zu optimized=%zu\n",
            tested,check->reasons[IC_FAN_PASS],check->reasons[IC_FAN_SUPPORT],
            check->reasons[IC_FAN_SEED],check->reasons[IC_FAN_SOLVE],check->reasons[IC_FAN_GUARD],
            check->reasons[IC_FAN_METRIC],check->reasons[IC_FAN_DIRECTION],
            check->min_sigma,check->max_sigma,check->max_rms,
            check->seed_certificates,check->optimized_certificates);
    if (check->reasons[IC_FAN_GUARD])
        fprintf(stderr,"[ink_closeup] single-slit guard refusals: reverted=%zu pure-SD=%zu total=%zu embedding=%zu\n",
                check->guard_reverted,check->guard_metric,check->guard_objective,check->guard_embedding);
    fprintf(stderr,"[ink_closeup] single-slit search: ");
    for (int k=0;k<CLOSEUP_CONE_FAN_WEIGHT_STAGES;k++)
        fprintf(stderr,"weight-%.9g=%zu ",ldexp(CLOSEUP_FIBER_DIRECTION_WEIGHT,-k),check->weight_certificates[k]);
    fprintf(stderr,"metric-refused-with-direction-pass=%zu\n",check->metric_refused_with_direction_pass);
}

#include "ink_closeup_cone_feasibility_test.inc"
