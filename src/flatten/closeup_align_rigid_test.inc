static int car_reject(void *context,size_t nv,const double *uv)
{ (void)nv;(void)uv;(*(int *)context)++;return 0; }

/* The separate rigid pass still requires a solver callback contract, but
 * must never call it: zero-iteration native guards certify the rigid map. */
static int car_unused(Arena_T a,void *c,size_t nv,size_t ne,const double *d,
    const int32_t *ea,const int32_t *eb,const double *blocks,const double *rhs,double *out)
{
    (void)a;(void)nv;(void)ne;(void)d;(void)ea;(void)eb;(void)blocks;(void)rhs;(void)out;
    (*(int *)c)++;return -1;
}

int CloseupAlign_rigid_selftest(void)
{
    const float xyz[12]={0,0,0,4,0,0,4,3,0,0,3,0};
    const int32_t faces[6]={0,1,2,0,2,3};
    const double original[8]={0,0,4,0,4,3,0,3};
    SlimRefineDirection arms[2]={{0},{0}};
    Arena_T arena=Arena_new_sized(4096);int failures=0,linear_calls=0,cases=0;
    CloseupAlignOpts opts;CloseupAlign_defaults(&opts);
    opts.slim.linear_solve=car_unused;opts.slim.linear_context=&linear_calls;
    opts.slim.verbose=0;opts.slim.threads=1;
    double uv[8],out[8];CloseupAlignStats report;
    for (int reflected=0;reflected<2;reflected++) for (int shifted=0;shifted<2;shifted++)
    for (int arm_variant=0;arm_variant<4;arm_variant++) {
        double angle=.61,co=cos(angle),si=sin(angle);
        for (int f=0;f<2;f++) {
            memset(arms+f,0,sizeof *arms);arms[f].weight=1;
            arms[f].h[0]=co;arms[f].h[1]=si;arms[f].v[0]=-si;arms[f].v[1]=co;
            if (arm_variant&1) for (int k=0;k<3;k++) arms[f].h[k]=-arms[f].h[k];
            if (arm_variant&2) for (int k=0;k<3;k++) {double t=arms[f].h[k];arms[f].h[k]=arms[f].v[k];arms[f].v[k]=t;}
        }
        for (int v=0;v<4;v++) {
            uv[2*v]=original[2*v]+(shifted?1e6:0);
            uv[2*v+1]=(reflected?-1:1)*original[2*v+1]+(shifted?-2e6:0);
        }
        int rc=CloseupAlign_rotate(arena,xyz,4,faces,2,uv,arms,NULL,&opts,out,&report);
        int ok=rc==CLOSEUP_ALIGN_PROPOSED && report.selected_seed==2 && report.target_met &&
            report.after.fiber_rms_degrees<1e-6 && fabs(report.after.sd_mean-report.before.sd_mean)<1e-8;
        for (int a=0;a<4;a++) for (int b=0;b<a;b++) {
            double old=hypot(uv[2*a]-uv[2*b],uv[2*a+1]-uv[2*b+1]);
            double now=hypot(out[2*a]-out[2*b],out[2*a+1]-out[2*b+1]);
            ok=ok && fabs(old-now)<1e-8;
        }
        if (!ok) failures++;cases++;
    }
    int rejected=0;opts.candidate_context=&rejected;opts.candidate_guard=car_reject;
    int rc=CloseupAlign_rotate(arena,xyz,4,faces,2,uv,arms,NULL,&opts,out,&report);
    if (rc!=CLOSEUP_ALIGN_UNCHANGED || !rejected || memcmp(out,uv,sizeof uv)) failures++;
    opts.candidate_context=NULL;opts.candidate_guard=NULL;
    arms[0].weight=arms[1].weight=0;
    rc=CloseupAlign_rotate(arena,xyz,4,faces,2,uv,arms,NULL,&opts,out,&report);
    if (rc!=CLOSEUP_ALIGN_NO_EVIDENCE || memcmp(out,uv,sizeof uv) || isfinite(report.after.fiber_rms_degrees)) failures++;
    /* The sweep's optimum is independently checked against the public angular
     * metric on a dense quarter-turn grid, including branch crossings. */
    uint32_t random=0x2934c913u;
    for (int test=0;test<32;test++) {
        Arena_Mark mark=Arena_save(arena);CaDomain d={arena,xyz,faces,arms,NULL,4,2,NULL};
        d.face=ARENA_CALLOC(arena,2,sizeof *d.face);
        for (int f=0;f<2;f++) {
            memset(arms+f,0,sizeof *arms);
            for (int a=0;a<2;a++) {
                random=1664525u*random+1013904223u;
                double theta=6.283185307179586*(double)(random&65535)/65536;
                double *arm=a?arms[f].v:arms[f].h;arm[0]=cos(theta);arm[1]=sin(theta);
            }
            arms[f].weight=.01+(f+1)*(test+1)*.13;
        }
        int ok=!ca_setup(&d,original);double theta=0;ok=ok && !ca_rigid_angle(&d,original,&theta);
        double co=cos(theta),si=sin(theta);CloseupAlignMetric best={0};
        for (int v=0;v<4;v++) {out[2*v]=co*original[2*v]-si*original[2*v+1];out[2*v+1]=si*original[2*v]+co*original[2*v+1];}
        ok=ok && !ca_measure(&d,out,1,&best);
        for (int sample=0;ok && sample<=720;sample++) {
            double angle=-.7853981633974483+1.5707963267948966*sample/720;
            co=cos(angle);si=sin(angle);CloseupAlignMetric m={0};
            for (int v=0;v<4;v++) {out[2*v]=co*original[2*v]-si*original[2*v+1];out[2*v+1]=si*original[2*v]+co*original[2*v+1];}
            ok=!ca_measure(&d,out,1,&m) && best.fiber_rms_degrees<=m.fiber_rms_degrees+1e-9;
        }
        if (!ok) failures++;cases++;Arena_restore(arena,mark);
    }
    if (linear_calls) failures++;
    Arena_dispose(&arena);
    fprintf(stderr,"[rigid selftest] %d controls, %d failures, %d nonlinear solves\n",cases+2,failures,linear_calls);
    return failures;
}
