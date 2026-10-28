typedef struct PMALossTest {double angle;int calls,force_cap,fail_at;} PMALossTest;

static int pmalt_linear(Arena_T arena,void *context,size_t nv,size_t ne,const double *diagonal,
    const int32_t *edge_a,const int32_t *edge_b,const double *blocks,const double *rhs,double *out)
{
    (void)arena;(void)ne;(void)diagonal;(void)edge_a;(void)edge_b;(void)blocks;(void)rhs;
    PMALossTest *t=context;t->calls++;if(nv!=4||(t->fail_at&&t->calls==t->fail_at))return -1;
    double error=atan2(out[3]-out[1],out[2]-out[0])+t->angle;
    double c=cos(-.5*error),s=sin(-.5*error);
    for(size_t v=0;v<nv;v++){double x=out[2*v],y=out[2*v+1];out[2*v]=c*x-s*y;out[2*v+1]=s*x+c*y;}
    return 0;
}
static int pmalt_material(void *context,const double *uv,size_t *pairs)
{(void)context;(void)uv;*pairs=0;return 1;}
static int pmalt_observe(void *context,const double *uv,const PMAEvent *e,const PMAMetric *m,size_t groups)
{(void)context;(void)uv;(void)e;(void)m;(void)groups;return 0;}
static int pmalt_refine(void *context,const double *current,const double *trial,const PMAEvent *e,const PMAMetric *m,size_t groups)
{(void)current;(void)trial;(void)e;(void)m;(void)groups;return ((PMALossTest *)context)->force_cap;}

int CloseupAngularLossStop_selftest(void)
{
    int failures=0,cases=0;PMAOptions o={0};o.loss_window=5;o.loss_relative_tolerance=1e-6;
    for(int sample=0;sample<32;sample++) {
        double area=1+10*sample,loss=100+.7*sample;PMAStats s={0};
        failures+=pmal_begin(&s,4*area+loss,area);cases++;
        /* A large SD constant cannot manufacture a small relative loss. */
        double next=loss*(1-1e-4);
        failures+=pmal_update(&s,&o,4*area+next,-loss,1)!=0||s.loss_stable_steps;cases++;
        /* Tiny accepted gain with a strong modeled descent is a stall. */
        for(int k=0;k<6;k++){next-=1e-9;failures+=pmal_update(&s,&o,4*area+next,-1,1)!=0||s.loss_stable_steps;cases++;}
        /* An unresolved inner cap cannot certify convergence either. */
        for(int k=0;k<6;k++){next-=1e-9;failures+=pmal_update(&s,&o,4*area+next,-1e-9,0)!=0||s.loss_stable_steps;cases++;}
        for(int k=0;k<5;k++) {
            next-=1e-9;int rc=pmal_update(&s,&o,4*area+next,-1e-9,1);
            failures+=rc!=(k==4)||s.loss_stable_steps!=k+1;cases++;
        }
        failures+=!s.loss_converged||s.termination!=PMA_STOP_LOSS_TOLERANCE;cases++;
        PMAStats before=s;
        failures+=pmal_update(&s,&o,NAN,0,1)!=-1||memcmp(&before,&s,sizeof s);cases++;
        failures+=pmal_update(&s,&o,4*area+next+1,0,1)!=-1||memcmp(&before,&s,sizeof s);cases++;
        failures+=pmal_update(&s,&o,4*area+next,NAN,1)!=-1||memcmp(&before,&s,sizeof s);cases++;
        failures+=pmal_begin(&s,4*area-1,area)!=-1;cases++;
    }
    /* Long-run policy: a short stretch of small gains must not stop.
     * Large model slopes, larger real gains and incomplete inner work each
     * prevent the 100-step certificate independently. */
    for(int scenario=0;scenario<4;scenario++) {
        PMAOptions long_options={0};PMAStats long_stats={0};
        long_options.loss_window=100;long_options.loss_relative_tolerance=1e-8;
        double value=112;
        failures+=pmal_begin(&long_stats,value,3);cases++;
        for(int k=0;k<100;k++) {
            value-=scenario==2?5e-6:1e-7;
            int rc=pmal_update(&long_stats,&long_options,value,scenario==1?-1e-3:-1e-9,scenario!=3);
            failures+=rc!=(scenario==0&&k==99)||long_stats.loss_stable_steps!=(scenario==0?k+1:0);cases++;
        }
        failures+=long_stats.loss_converged!=(scenario==0);cases++;
    }
    const float xyz[12]={0,0,0,1,0,0,1,1,0,0,1,0};const int32_t faces[6]={0,1,2,0,2,3};
    const double seed[8]={0,0,1,0,1,1,0,1};uint8_t mask[2]={1,1},guard[1]={0};
    for(int sample=0;sample<16;sample++)for(int scenario=0;scenario<4;scenario++) {
        Arena_T arena=Arena_new();PMALossTest t={0};t.angle=.08+.01*sample;
        t.force_cap=scenario==2;t.fail_at=scenario==3?1:0;
        SlimRefineDirection arms[2]={0};double out[8];PMAStats stats={0};
        for(int f=0;f<2;f++){arms[f].h[0]=cos(t.angle);arms[f].h[1]=sin(t.angle);arms[f].v[0]=-sin(t.angle);arms[f].v[1]=cos(t.angle);arms[f].weight=1;}
        PMAOptions opts={0};SlimRefine_defaults(&opts.slim);opts.slim.strip_pack=0;opts.slim.guard_boundary=1;
        opts.slim.directions=arms;opts.slim.direction_weight=16;opts.slim.linear_context=&t;opts.slim.linear_solve=pmalt_linear;
        opts.model=1;opts.iterations=scenario==1?2:80;opts.refinement_limit=1;
        opts.groups=1;opts.group_masks=mask;opts.group_count_guard=guard;opts.reference_uv=seed;opts.context=&t;
        opts.material_guard=pmalt_material;opts.observe=pmalt_observe;opts.refine=pmalt_refine;
        opts.loss_window=5;opts.loss_relative_tolerance=1e-6;
        int rc=CloseupAngular_run(arena,xyz,4,faces,2,seed,&opts,out,&stats),bad=0;
        if(scenario==0)bad=rc||!stats.loss_converged||stats.iterations>=80||stats.loss_stable_steps!=5||
            stats.termination!=PMA_STOP_LOSS_TOLERANCE||stats.last_relative_loss_decrease>1e-6||stats.last_relative_model_slope>1e-6;
        if(scenario==1)bad=rc||stats.loss_converged||stats.iterations!=2||stats.termination!=PMA_STOP_WORK_LIMIT;
        if(scenario==2)bad=rc||stats.loss_converged||!stats.refinement_limit_hits||stats.loss_stable_steps;
        if(scenario==3)bad=rc==0||stats.loss_converged||stats.termination!=PMA_STOP_SOLVER_FAILURE||memcmp(out,seed,sizeof out);
        if(bad)fprintf(stderr,"[loss control FAIL] sample=%d scenario=%d rc=%d finish=%d iterations=%d stable=%d loss=%.9g slope=%.9g\n",sample,scenario,rc,stats.termination,stats.iterations,stats.loss_stable_steps,stats.final_loss,stats.last_relative_model_slope);
        failures+=bad;cases++;Arena_dispose(&arena);
    }
    fprintf(stderr,"[loss stopping selftest] %d cases, %d failures\n",cases,failures);return failures?1:0;
}
