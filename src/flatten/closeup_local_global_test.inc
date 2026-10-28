typedef struct ClgTest { int calls,refuse; } ClgTest;
static int clg_test_guard(void *context,const double *uv,size_t *pairs)
{
    ClgTest *s=context;(void)uv;*pairs=0;s->calls++;
    return s->refuse?0:1;
}

static int clg_test(void)
{
    const float xyz[12]={0,0,0,1,0,0,1,1,0,0,1,0};
    const int32_t faces[6]={0,1,2,0,2,3};
    const double reference[8]={0,0,1,0,1,1,0,1};
    const uint8_t masks[2]={1,1},flags[1]={1};
    int failed=0,cases=0;
    for(int test=0;test<4;test++) {
        Arena_T arena=Arena_new();PMAOptions o={0};PMAStats stats={0};ClgTest guard={0};
        SlimRefine_defaults(&o.slim);SlimRefineDirection arms[2]={0};
        double angle=test==3?0:.35,c=cos(angle),s=sin(angle),out[8],seed[8];
        memcpy(seed,reference,sizeof seed);
        for(size_t f=0;f<2;f++) {
            arms[f].h[0]=c;arms[f].h[1]=s;arms[f].v[0]=-s;arms[f].v[1]=c;arms[f].weight=1;
        }
        o.slim.directions=arms;o.slim.direction_weight=16;o.slim.direction_rigid_gauge=0;
        o.slim.strip_pack=0;o.slim.guard_boundary=1;o.iterations=test==1?1:50;
        o.refinement_limit=64;o.loss_window=100;o.loss_relative_tolerance=1e-8;
        o.target_rms_degrees=test==1?.00001:10;
        o.reference_uv=reference;o.groups=1;o.group_masks=masks;o.group_count_guard=flags;
        o.context=&guard;o.material_guard=clg_test_guard;guard.refuse=test==2;
        int status=CloseupLocalGlobal_run(arena,xyz,4,faces,2,seed,&o,out,&stats);
        int ok=!memcmp(seed,reference,sizeof seed);
        if(test==2)ok&=status<0&&!memcmp(out,seed,sizeof out)&&stats.termination==PMA_STOP_SOLVER_FAILURE;
        else if(test==1)ok&=!status&&!stats.loss_converged&&stats.termination==PMA_STOP_WORK_LIMIT&&stats.accepted_steps==1;
        else if(test==3)ok&=!status&&!stats.accepted_steps&&stats.termination==PMA_STOP_TARGET&&!memcmp(out,seed,sizeof out);
        else ok&=!status&&stats.accepted_steps>0&&stats.final_loss<stats.initial_loss&&stats.termination==PMA_STOP_TARGET;
        fprintf(stderr,"[closeup local global test] case=%d status=%d termination=%d accepted=%d %s\n",
            test,status,stats.termination,stats.accepted_steps,ok?"PASS":"FAIL");
        failed+=!ok;cases++;Arena_dispose(&arena);
    }
    fprintf(stderr,"closeup local global adapter: %d cases, %d failures\n",cases,failed);
    return failed?1:0;
}
