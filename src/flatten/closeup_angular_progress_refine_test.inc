int CloseupAngularProgressRefine_selftest(void)
{
    int failures=0,cases=0;
    for(int sample=0;sample<32;sample++) {
        PMAOptions o={0};o.loss_window=5;o.loss_relative_tolerance=1e-6;
        double area=1+sample,objective=4*area+10,derivative=-1;
        failures+=pma_progress_restart(&o,1,objective,area,derivative)!=0;cases++;
        o.refinement_progress_alpha=.5;
        failures+=pma_progress_restart(&o,.25,objective,area,derivative)!=0;cases++;
        failures+=pma_progress_restart(&o,.5,objective,area,derivative)!=1;cases++;
        failures+=pma_progress_restart(&o,1,objective,area,derivative)!=1;cases++;
        failures+=pma_progress_restart(&o,1,objective,area,-1e-8)!=0;cases++;
        failures+=pma_progress_restart(&o,1,objective,area,NAN)!=0;cases++;
        failures+=pma_progress_restart(&o,NAN,objective,area,derivative)!=0;cases++;
        failures+=pma_progress_restart(&o,1,INFINITY,area,derivative)!=0;cases++;
    }
    const float xyz[12]={0,0,0,1,0,0,1,1,0,0,1,0};const int32_t faces[6]={0,1,2,0,2,3};
    const double seed[8]={0,0,1,0,1,1,0,1};uint8_t mask[2]={1,1},guard[1]={0};
    for(int sample=0;sample<16;sample++)for(int scenario=0;scenario<4;scenario++) {
        Arena_T arena=Arena_new();PMALossTest t={0};t.angle=.08+.01*sample;
        t.force_cap=scenario!=1;t.fail_at=scenario==3?1:0;
        SlimRefineDirection arms[2]={0};double out[8];PMAStats stats={0};
        for(int f=0;f<2;f++){arms[f].h[0]=cos(t.angle);arms[f].h[1]=sin(t.angle);arms[f].v[0]=-sin(t.angle);arms[f].v[1]=cos(t.angle);arms[f].weight=1;}
        PMAOptions o={0};SlimRefine_defaults(&o.slim);o.slim.strip_pack=0;o.slim.guard_boundary=1;
        o.slim.directions=arms;o.slim.direction_weight=16;o.slim.linear_context=&t;o.slim.linear_solve=pmalt_linear;
        o.model=1;o.iterations=scenario==0?2:80;o.refinement_limit=64;o.refinement_progress_alpha=.5;
        o.groups=1;o.group_masks=mask;o.group_count_guard=guard;o.reference_uv=seed;o.context=&t;
        o.material_guard=pmalt_material;o.observe=pmalt_observe;o.refine=pmalt_refine;
        o.loss_window=5;o.loss_relative_tolerance=1e-6;
        int rc=CloseupAngular_run(arena,xyz,4,faces,2,seed,&o,out,&stats),bad=0;
        if(scenario==0)bad=rc||stats.accepted_steps!=2||t.calls!=2||stats.refinement_progress_stops!=2||stats.refinement_limit_hits||stats.loss_converged||stats.loss_stable_steps;
        if(scenario==1)bad=rc||!stats.loss_converged||stats.refinement_progress_stops||stats.loss_stable_steps!=5;
        if(scenario==2)bad=rc||stats.loss_converged||!stats.refinement_progress_stops||!stats.refinement_limit_hits||stats.loss_stable_steps;
        if(scenario==3)bad=!rc||stats.loss_converged||memcmp(out,seed,sizeof out);
        if(bad)fprintf(stderr,"[progress control FAIL] sample=%d scenario=%d rc=%d steps=%d progress=%zu caps=%zu stable=%d\n",sample,scenario,rc,stats.accepted_steps,stats.refinement_progress_stops,stats.refinement_limit_hits,stats.loss_stable_steps);
        failures+=bad;cases++;Arena_dispose(&arena);
    }
    fprintf(stderr,"[progress refinement selftest] %d cases, %d failures\n",cases,failures);return failures?1:0;
}
