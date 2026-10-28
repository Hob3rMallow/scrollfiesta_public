int CloseupAngularGlobalRefine_selftest(void)
{
    int failures=0,cases=0;const int limits[3]={4,16,64};
    const float xyz[12]={0,0,0,1,0,0,1,1,0,0,1,0};const int32_t faces[6]={0,1,2,0,2,3};
    const double seed[8]={0,0,1,0,1,1,0,1};uint8_t mask[2]={1,1},guard[1]={0};
    for(int sample=0;sample<16;sample++)for(int steps=1;steps<=3;steps++) {
        double baseline[8]={0};
        for(int variant=0;variant<3;variant++) {
            Arena_T arena=Arena_new();PMALossTest t={0};t.angle=.08+.01*sample;t.force_cap=1;
            SlimRefineDirection arms[2]={0};double out[8];PMAStats stats={0};
            for(int f=0;f<2;f++){arms[f].h[0]=cos(t.angle);arms[f].h[1]=sin(t.angle);arms[f].v[0]=-sin(t.angle);arms[f].v[1]=cos(t.angle);arms[f].weight=1;}
            PMAOptions o={0};SlimRefine_defaults(&o.slim);o.slim.strip_pack=0;o.slim.guard_boundary=1;
            o.slim.directions=arms;o.slim.direction_weight=16;o.slim.linear_context=&t;o.slim.linear_solve=pmalt_linear;
            o.model=1;o.iterations=steps;o.refinement_limit=limits[variant];
            o.groups=1;o.group_masks=mask;o.group_count_guard=guard;o.reference_uv=seed;o.context=&t;
            o.material_guard=pmalt_material;o.observe=pmalt_observe;o.refine=pmalt_refine;
            int rc=CloseupAngular_run(arena,xyz,4,faces,2,seed,&o,out,&stats),expected=steps*(limits[variant]+1);
            if(!variant)memcpy(baseline,out,sizeof baseline);
            failures+=rc||t.calls!=expected||stats.accepted_steps!=steps||stats.refinement_solves!=(size_t)expected||
                stats.refinement_limit_hits!=(size_t)steps||memcmp(baseline,out,sizeof baseline);
            cases++;Arena_dispose(&arena);
        }
    }
    fprintf(stderr,"[global refinement selftest] %d cases, %d failures\n",cases,failures);return failures?1:0;
}
