typedef struct PMARefineTest {
    int scenario,calls,refines,observations,commits,material_calls,matrix_mismatch;
    double target[3][8],last[8],matrix[64];size_t matrix_size;
} PMARefineTest;

static int pmart_linear(Arena_T arena,void *context,size_t nv,size_t ne,const double *diagonal,
    const int32_t *edge_a,const int32_t *edge_b,const double *blocks,const double *rhs,double *out)
{
    (void)arena;(void)edge_a;(void)edge_b;PMARefineTest *t=context;int call=t->calls++;
    if(nv!=4||3*nv+3*ne+2*nv>64)return -1;
    double matrix[64];size_t n=3*nv+3*ne+2*nv;
    memcpy(matrix,diagonal,3*nv*sizeof *matrix);memcpy(matrix+3*nv,blocks,3*ne*sizeof *matrix);
    memcpy(matrix+3*nv+3*ne,rhs,2*nv*sizeof *matrix);
    if(!call){memcpy(t->matrix,matrix,n*sizeof *matrix);t->matrix_size=n;}
    else t->matrix_mismatch|=n!=t->matrix_size||memcmp(t->matrix,matrix,n*sizeof *matrix)!=0;
    if(t->scenario==10||(t->scenario==4&&call))return -1;
    memcpy(out,t->target[call<3?call:2],8*sizeof *out);
    if(t->scenario==6&&call)out[0]=NAN;
    return 0;
}
static int pmart_material(void *context,const double *uv,size_t *pairs)
{
    (void)uv;PMARefineTest *t=context;t->material_calls++;*pairs=0;
    if(t->material_calls==4&&t->scenario==7)return -1;
    if(t->material_calls==4&&t->scenario==12){*pairs=1;return 0;}
    return 1;
}
static int pmart_refine(void *context,const double *current,const double *trial,
    const PMAEvent *e,const PMAMetric *m,size_t groups)
{
    (void)current;(void)trial;(void)e;(void)m;(void)groups;PMARefineTest *t=context;t->refines++;
    if(t->scenario==8)return -1;
    if(t->scenario==0||t->scenario==11)return 0;
    if(t->scenario==3)return 1;
    return t->refines<(t->scenario==2?3:2);
}
static int pmart_observe(void *context,const double *uv,const PMAEvent *e,const PMAMetric *m,size_t groups)
{
    (void)m;(void)groups;PMARefineTest *t=context;t->observations++;
    if((t->scenario==9&&e->deferred)||(t->scenario==11&&e->accepted))return 1;
    if(e->accepted){t->commits++;memcpy(t->last,uv,8*sizeof *uv);}
    return 0;
}

int CloseupAngularRefineStep_selftest(void)
{
    int failures=0,cases=0;
    const float xyz[12]={0,0,0,1,0,0,1,1,0,0,1,0};const int32_t faces[6]={0,1,2,0,2,3};
    const double seed[8]={0,0,1,0,1,1,0,1};uint8_t mask[2]={1,1},guard[1]={0};
    for(int sample=0;sample<16;sample++)for(int scenario=0;scenario<13;scenario++) {
        PMARefineTest t={0};t.scenario=scenario;double angle=.15+.01*sample,out[8];SlimRefineDirection arms[2]={0};
        for(int f=0;f<2;f++){arms[f].h[0]=cos(angle);arms[f].h[1]=sin(angle);arms[f].v[0]=-sin(angle);arms[f].v[1]=cos(angle);arms[f].weight=1;}
        for(int a=0;a<3;a++) {
            double theta=-angle*(a==0?1.0/3:a==1?2.0/3:scenario==3?.8:1.0/6);
            for(int v=0;v<4;v++){t.target[a][2*v]=cos(theta)*seed[2*v]-sin(theta)*seed[2*v+1];t.target[a][2*v+1]=sin(theta)*seed[2*v]+cos(theta)*seed[2*v+1];}
        }
        if(scenario==5)memcpy(t.target[1],seed,sizeof seed);
        PMAOptions o={0};SlimRefine_defaults(&o.slim);o.slim.strip_pack=0;o.slim.guard_boundary=1;
        o.slim.directions=arms;o.slim.direction_weight=16;o.slim.linear_solve=pmart_linear;o.slim.linear_context=&t;
        o.model=PMA_JOINT_ANGLE_V1;o.iterations=1;o.groups=1;o.group_masks=mask;o.group_count_guard=guard;o.reference_uv=seed;
        o.context=&t;o.material_guard=pmart_material;o.observe=pmart_observe;o.refine=pmart_refine;o.refinement_limit=2;
        Arena_T arena=Arena_new();PMAStats stats={0};int rc=CloseupAngular_run(arena,xyz,4,faces,2,seed,&o,out,&stats);
        int refused=scenario==7||scenario==9||scenario==10||scenario==11||scenario==12;
        int error=scenario==9||scenario==10||scenario==11;
        int best=scenario==1||scenario==2?1:scenario==3?2:0;
        int bad=(error?rc>=0:rc!=0)||t.matrix_mismatch;
        if(refused)bad|=memcmp(out,seed,sizeof out)!=0||stats.accepted_steps!=0||t.commits!=0;
        else bad|=memcmp(out,t.target[best],sizeof out)!=0||memcmp(out,t.last,sizeof out)!=0||stats.accepted_steps!=1||t.commits!=1;
        if(scenario==2)bad|=stats.refinement_prior_best!=1||t.calls!=3;
        if(scenario==3)bad|=stats.refinement_limit_hits!=1||t.calls!=3;
        if(scenario==4||scenario==5||scenario==6)bad|=stats.refinement_failures!=1||stats.refinement_prior_best!=1||t.calls!=2;
        if(scenario==8)bad|=stats.refinement_failures!=1||t.calls!=1;
        if(scenario==7||scenario==12)bad|=stats.stop_reason!=(scenario==7?PMA_MATERIAL_UNKNOWN:PMA_MATERIAL);
        if(bad)fprintf(stderr,"[same quadratic control FAIL] sample=%d scenario=%d rc=%d accepted=%d calls=%d prior=%zu failures=%zu stop=%d\n",sample,scenario,rc,stats.accepted_steps,t.calls,stats.refinement_prior_best,stats.refinement_failures,stats.stop_reason);
        failures+=bad;cases++;Arena_dispose(&arena);
    }
    fprintf(stderr,"[same quadratic controls] %d cases, %d failures\n",cases,failures);return failures?1:0;
}
