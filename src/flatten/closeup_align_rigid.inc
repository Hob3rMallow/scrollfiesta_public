typedef struct CaRigidEvent {
    double crossing, mean, weight;
    size_t face;
} CaRigidEvent;

static int ca_rigid_compare(const void *a,const void *b)
{
    const CaRigidEvent *x=a,*y=b;
    if (x->crossing!=y->crossing) return x->crossing<y->crossing?-1:1;
    return x->face<y->face?-1:x->face>y->face;
}

static void ca_rigid_interval(double mass,double linear,double square,
                                double lo,double hi,double *best,double *angle)
{
    double t=fmax(lo,fmin(hi,-linear/mass));
    double value=square+2*t*linear+t*t*mass;
    if (value<*best) { *best=value;*angle=t; }
}

/* Exact joint angular score, not the mapped-length directional objective.
 * With mean=(angle(h)+angle(v)-pi/2)/2 and delta=angle(h)-angle(v)+pi/2,
 * the squared joint error is wrap_pi/2(mean+turn)^2+wrap_pi(delta)^2/4.
 * The second term is fixed. Sorted branch crossings make the first term a
 * finite sequence of scalar quadratics. This proposes only one rigid turn. */
static int ca_rigid_angle(CaDomain *d,const double *uv,double *angle)
{
    const double period=1.57079632679489661923;
    Arena_Mark mark=Arena_save(d->arena);
    CaRigidEvent *events=ARENA_ALLOC(d->arena,d->nf*sizeof *events);
    size_t count=0;double mass=0,linear=0,square=0;
    for (size_t f=0;f<d->nf;f++) {
        const CaFace *q=d->face+f;
        if (!(q->evidence>0)) continue;
        double j[4];ca_jacobian(d,f,uv,j);
        const double *a=q->frame;
        double hx=j[0]*a[0]+j[1]*a[2],hy=j[2]*a[0]+j[3]*a[2];
        double vx=j[0]*a[1]+j[1]*a[3],vy=j[2]*a[1]+j[3]*a[3];
        double mean=.5*(atan2(hy,hx)+atan2(vy,vx)-period);
        mean-=period*floor((mean+.5*period)/period);
        double w=q->area*q->evidence;
        events[count].mean=mean;events[count].crossing=.5*period-mean;
        events[count].weight=w;events[count].face=f;count++;
        mass+=w;linear+=w*mean;square+=w*mean*mean;
    }
    if (!(mass>0) || !isfinite(mass+linear+square)) { Arena_restore(d->arena,mark);return -1; }
    qsort(events,count,sizeof *events,ca_rigid_compare);
    double best=square,best_angle=0,lo=0;
    for (size_t e=0;e<count;e++) {
        const CaRigidEvent *v=events+e;
        ca_rigid_interval(mass,linear,square,lo,v->crossing,&best,&best_angle);
        linear-=v->weight*period;
        square+=v->weight*(period*period-2*period*v->mean);
        lo=v->crossing;
    }
    ca_rigid_interval(mass,linear,square,lo,period,&best,&best_angle);
    if (best_angle>.5*period) best_angle-=period;
    *angle=best_angle;
    Arena_restore(d->arena,mark);return isfinite(*angle)?0:-1;
}

static int ca_try_rigid(CaDomain *d,const CloseupAlignOpts *opts,const double *seed,
                          int seed_id,const double center[2],int orientation,
                          double *trial,double *best,CloseupAlignStats *report)
{
    double angle=0,origin[2];
    if (ca_rigid_angle(d,seed,&angle) || fabs(angle)<=1e-12) return 0;
    ca_centroid(d,seed,origin);
    /* A rejected full turn may still have a useful safe fraction. Each exact
     * endpoint is independently certified; contact feasibility is not assumed
     * monotone, and no rejected candidate changes the retained winner. */
    for (int attempt=0;attempt<CLOSEUP_ALIGN_RIGID_ATTEMPTS;attempt++,angle*=.5) {
        if (report->trials>=CLOSEUP_ALIGN_MAX_TRIALS) return -1;
        CloseupAlignTrial *record=report->trial+report->trials++;
        record->seed=seed_id;record->direction_weight=0;record->solver_status=0;
        record->candidate_guard_pass=opts->candidate_guard?0:1;
        record->metric.fiber_rms_degrees=NAN;
        double co=cos(angle),si=sin(angle);
        for (size_t v=0;v<d->nv;v++) {
            double x=seed[2*v]-origin[0],y=seed[2*v+1]-origin[1];
            trial[2*v]=origin[0]+co*x-si*y;
            trial[2*v+1]=origin[1]+si*x+co*y;
        }
        int rc=ca_output_candidate(d,opts,trial,center,orientation,best,report,record);
        if (opts->slim.verbose) fprintf(stderr,
            "[rigid trial] seed=%d attempt=%d degrees=%.12g feasible=%d metric=%d contact=%d selected=%d RMS=%.12g\n",
            seed_id,attempt,angle*57.29577951308232,rc,record->scale_nonregression,
            record->candidate_guard_pass,record->selected,record->metric.fiber_rms_degrees);
        if (record->selected) return 1;
    }
    return 0;
}
