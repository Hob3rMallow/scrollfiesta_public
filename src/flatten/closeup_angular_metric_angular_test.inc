#include "sparse_solve.h"
static int pma_test_solve(Arena_T arena,void *context,size_t nv,size_t ne,const double *diagonal,
    const int32_t *edge_a,const int32_t *edge_b,const double *blocks,const double *rhs,double *result)
{
    (void)context;Arena_Mark mark=Arena_save(arena);volatile int rc=-1;
    TRY {
        size_t count=3*nv+4*ne,at=0;int *row=ARENA_ALLOC(arena,count*sizeof *row),*col=ARENA_ALLOC(arena,count*sizeof *col);
        double *value=ARENA_ALLOC(arena,count*sizeof *value);
        for(size_t v=0;v<nv;v++) {
            int x=(int)(2*v),y=x+1;row[at]=col[at]=x;value[at++]=diagonal[3*v];
            row[at]=y;col[at]=x;value[at++]=diagonal[3*v+1];row[at]=col[at]=y;value[at++]=diagonal[3*v+2];
        }
        for(size_t e=0;e<ne;e++)for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
            int x=2*edge_a[e]+a,y=2*edge_b[e]+b;row[at]=x>y?x:y;col[at]=x<y?x:y;value[at++]=blocks[3*e+(size_t)(a+b)];
        }
        rc=Sparse_solve_sym((int)(2*nv),(int)count,row,col,value,rhs,result,SPARSE_SPD);
    } FINALLY {Arena_restore(arena,mark);} END_TRY;
    return rc;
}
typedef struct PMATestContext {size_t nv;int calls,accepted,physical_refusals,unknown;double *last,*target;double alpha;} PMATestContext;
static int pma_test_material(void *context,const double *uv,size_t *pairs)
{(void)uv;PMATestContext *c=context;*pairs=0;return c->unknown?-1:1;}
static int pma_test_observe(void *context,const double *uv,const PMAEvent *e,const PMAMetric *m,size_t groups)
{
    (void)m;(void)groups;PMATestContext *c=context;c->physical_refusals+=e->first_failure==PMA_PHYSICAL;
    if(e->accepted){c->accepted++;c->alpha=e->alpha;if(c->last)memcpy(c->last,uv,2*c->nv*sizeof *uv);}return 0;
}
static int pma_test_first_target(Arena_T arena,void *context,size_t nv,size_t ne,const double *d,
    const int32_t *a,const int32_t *b,const double *e,const double *rhs,double *out)
{
    (void)arena;(void)ne;(void)d;(void)a;(void)b;(void)e;(void)rhs;PMATestContext *c=context;
    if(!c->calls++)memcpy(out,c->target,2*nv*sizeof *out);return 0;
}
static void pma_test_options(PMAOptions *o,PMATestContext *context,const SlimRefineDirection *arms,
    const double *reference,const uint8_t *mask,const uint8_t *count_guard)
{
    memset(o,0,sizeof *o);SlimRefine_defaults(&o->slim);o->slim.strip_pack=0;o->slim.guard_boundary=1;
    o->slim.directions=arms;o->slim.direction_weight=16;o->slim.linear_solve=pma_test_solve;
    o->iterations=10;o->model=PMA_JOINT_ANGLE_V1;o->groups=1;o->group_masks=mask;o->group_count_guard=count_guard;
    o->reference_uv=reference;o->context=context;o->material_guard=pma_test_material;o->observe=pma_test_observe;
}
static int pma_expect(int condition,const char *name,int *checks)
{(*checks)++;fprintf(stderr,"[metric angular control] %s %s\n",condition?"PASS":"FAIL",name);return !condition;}

int CloseupAngular_selftest(void)
{
    int failures=PrivateJointAngle_selftest()+CloseupAngularRefineStep_selftest(),checks=0;
    {
        PMAMetric before={0},after={0};PMAOptions options={0};uint8_t guard=0;size_t failed=0;
        double capacity=1;before.faces=after.faces=100;before.area=after.area=100;
        options.groups=1;options.group_count_guard=&guard;options.group_bad10_limits=&capacity;
        after.bad10_area=2;
        failures+=pma_expect(!pma_physical(&after,&before,&options,&failed),
            "stationary components' reserved area tightens the local physical budget",&checks);
        after.bad10_area=1;
        failures+=pma_expect(pma_physical(&after,&before,&options,&failed),"remaining area capacity is usable",&checks);
        capacity=0;after.bad10_area=.001;
        failures+=pma_expect(!pma_physical(&after,&before,&options,&failed),"zero budget preserves original good-face subsets",&checks);
        capacity=NAN;
        failures+=pma_expect(!pma_physical(&after,&before,&options,&failed),"invalid area cap fails closed",&checks);
    }
    const int32_t triangle[3]={0,1,2};SRDirection d={{-1,1,0},{-1,0,1},16};
    double a=.3,c=cos(a),s=sin(a),u[6]={0,0,c,s,-s,c},half[6];for(int k=0;k<6;k++)half[k]=.5*u[k];
    double angle0=0,angle1=0,e0=sr_direction_energy(&d,triangle,u,&angle0),e1=sr_direction_energy(&d,triangle,half,&angle1);
    failures+=pma_expect(fabs(e1-.25*e0)<1e-12&&fabs(angle1-angle0)<1e-12&&fabs(sra_energy(&d,triangle,half)-sra_energy(&d,triangle,u))<1e-12,
        "radial compression reduces mapped length by 75% with identical angle; angular objective invariant",&checks);
    failures+=pma_expect(fabs(sra_derivative(&d,triangle,u,half))<1e-12,"angular radial derivative is zero",&checks);
    {
        float xyz[12]={0,0,0,1,0,0,1,1,0,0,1,0};int32_t faces[6]={0,1,2,0,2,3};
        double seed[8]={0,0,1,0,1,1,0,1},target[8],last[8],out[8];uint8_t mask[2]={1,1},cg[1]={0};SlimRefineDirection arms[2]={0};
        for(int f=0;f<2;f++){arms[f].h[0]=c;arms[f].h[1]=s;arms[f].v[0]=-s;arms[f].v[1]=c;arms[f].weight=1;}
        for(int v=0;v<4;v++){target[2*v]=1.2*(c*seed[2*v]+s*seed[2*v+1]);target[2*v+1]=1.2*(-s*seed[2*v]+c*seed[2*v+1]);}
        PMATestContext context={0};context.nv=4;context.last=last;context.target=target;PMAOptions o;PMAStats stats;
        pma_test_options(&o,&context,arms,seed,mask,cg);o.slim.linear_solve=pma_test_first_target;o.slim.linear_context=&context;o.iterations=3;
        Arena_T arena=Arena_new();int rc=CloseupAngular_run(arena,xyz,4,faces,2,seed,&o,out,&stats);
        failures+=pma_expect(!rc&&context.physical_refusals>0&&context.accepted==1&&context.alpha>0&&context.alpha<1,
            "full proposal violates metric but a uniform fraction is retained",&checks);
        failures+=pma_expect(!rc&&stats.iterations==2&&memcmp(out,last,sizeof out)==0&&memcmp(out,seed,sizeof out)!=0,
            "later stationary refusal preserves exact accepted intermediate",&checks);
        context.unknown=1;context.calls=0;
        rc=CloseupAngular_run(arena,xyz,4,faces,2,seed,&o,out,&stats);
        failures+=pma_expect(rc<0&&memcmp(out,seed,sizeof out)==0,"unknown material query refuses without candidate adoption",&checks);
        Arena_dispose(&arena);
    }
    {
        /* Two components translate through each other; both endpoint maps
         * are individually valid. Continuous all-piece rim proof must fail. */
        float xyz[24]={0,0,0,1,0,0,1,1,0,0,1,0,4,0,0,5,0,0,5,1,0,4,1,0};
        int32_t faces[12]={0,1,2,0,2,3,4,5,6,4,6,7};double uv[16],target[16];
        for(int v=0;v<8;v++){uv[2*v]=xyz[3*v];uv[2*v+1]=xyz[3*v+1];target[2*v]=uv[2*v]+(v<4?4:-4);target[2*v+1]=uv[2*v+1];}
        Arena_T arena=Arena_new();SlimRefineOpts o;SlimRefine_defaults(&o);SRSystem system;
        int rc=sr_build_system(arena,xyz,8,faces,4,&o,1,&system);SREdgeBox boxes[8];
        failures+=pma_expect(!rc&&sr_boundary_intersections(&system,target,boxes)==0&&pma_rim_path(&system,uv,target,boxes)>0,
            "cross-component intermediate crossing refused despite endpoint component validity",&checks);
        for(int k=0;k<16;k++)target[k]=uv[k]+.125*(target[k]-uv[k]);
        failures+=pma_expect(!rc&&!pma_rim_path(&system,uv,target,boxes),"safe uniform prefix certified across all components",&checks);
        Arena_dispose(&arena);
    }
    {
        enum {NX=21,NY=5,NV=NX*NY,NF=2*(NX-1)*(NY-1)};
        float xyz[3*NV];int32_t faces[3*NF];double flat[2*NV],bent[2*NV],out[2*NV],last[2*NV];
        SlimRefineDirection arms[NF];uint8_t mask[NF],cg[1]={0};memset(arms,0,sizeof arms);memset(mask,1,sizeof mask);
        for(int y=0;y<NY;y++)for(int x=0;x<NX;x++) {
            int v=y*NX+x;double xx=-100+10*x,yy=-10+5*y,t=xx/200,r=200-yy;
            xyz[3*v]=(float)xx;xyz[3*v+1]=(float)yy;xyz[3*v+2]=0;flat[2*v]=xx;flat[2*v+1]=yy;
            bent[2*v]=r*sin(t);bent[2*v+1]=200-r*cos(t);
        }
        for(int y=0,f=0;y<NY-1;y++)for(int x=0;x<NX-1;x++) {
            int v=y*NX+x;faces[3*f]=v;faces[3*f+1]=v+1;faces[3*f+2]=v+NX+1;f++;
            faces[3*f]=v;faces[3*f+1]=v+NX+1;faces[3*f+2]=v+NX;f++;
        }
        for(int f=0;f<NF;f++) {
            int aa=faces[3*f],b=faces[3*f+1],cc=faces[3*f+2];
            double x1=flat[2*b]-flat[2*aa],y1=flat[2*b+1]-flat[2*aa+1],x2=flat[2*cc]-flat[2*aa],y2=flat[2*cc+1]-flat[2*aa+1];
            double v1=bent[2*b+1]-bent[2*aa+1],v2=bent[2*cc+1]-bent[2*aa+1],det=x1*y2-x2*y1;
            double j10=(v1*y2-v2*y1)/det,j11=(-v1*x2+v2*x1)/det,n=hypot(j10,j11);
            arms[f].h[0]=j11/n;arms[f].h[1]=-j10/n;arms[f].v[0]=-arms[f].h[1];arms[f].v[1]=arms[f].h[0];arms[f].weight=1;
        }
        Arena_T arena=Arena_new();PMATestContext context={0};context.nv=NV;context.last=last;PMAOptions o;PMAStats stats;
        pma_test_options(&o,&context,arms,flat,mask,cg);
        int rc=CloseupAngular_run(arena,xyz,NV,faces,NF,flat,&o,out,&stats);SRSystem system;
        PMAMetric before[1],after[1];int valid=!sr_build_system(arena,xyz,NV,faces,NF,&o.slim,1,&system)&&!pma_measure(&system,flat,&o,before)&&!pma_measure(&system,out,&o,after);
        size_t failed=0,moved=0;for(int v=0;v<NV;v++)if((v%NX==0||v%NX==NX-1||v<NX||v>=NV-NX)&&hypot(out[2*v]-flat[2*v],out[2*v+1]-flat[2*v+1])>1e-8)moved++;
        failures+=pma_expect(!rc&&valid&&stats.accepted_steps>0&&after[0].angular16<.5*before[0].angular16&&pma_physical(after,before,&o,&failed)&&moved==2*NX+2*NY-4,
            "nonconstant physical cross native TAUCS angular solve improves with every boundary vertex free",&checks);
        fprintf(stderr,"[metric angular strip] accepted=%d SD %.12g -> %.12g RMS %.12g -> %.12g boundary_moved=%zu/%d\n",
            stats.accepted_steps,before[0].sd,after[0].sd,sqrt(before[0].angular16/(32*before[0].mass))*180/3.14159265358979323846,
            sqrt(after[0].angular16/(32*after[0].mass))*180/3.14159265358979323846,moved,2*NX+2*NY-4);
        Arena_dispose(&arena);
    }
    fprintf(stderr,"[metric angular controls] checks=%d failures=%d\n",checks,failures);return failures?1:0;
}
