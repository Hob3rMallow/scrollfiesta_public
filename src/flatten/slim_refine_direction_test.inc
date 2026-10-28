/* Native tests for physical, unsigned material-frame constraints. */
static int sr_direction_selftest(void)
{
    int failures=0;
    {
        /* One common turn must preserve placement even when components have
         * opposite input handedness. Zero iterations isolates the gauge. */
        float xyz[24];const int32_t faces[12]={0,1,2,0,2,3,4,5,6,4,6,7};
        const double corners[8]={0,0,1,0,1,1,0,1};double seed[16],result[16],center[2]={0,0};
        SlimRefineDirection dirs[4];double theta=.37;
        for (size_t v=0;v<8;v++) {
            double x=corners[2*(v%4)]+(v<4?-3:5),y=(v<4?1:-1)*corners[2*(v%4)+1]+(v<4?0:2);
            xyz[3*v]=(float)x;xyz[3*v+1]=(float)corners[2*(v%4)+1];xyz[3*v+2]=0;
            seed[2*v]=cos(theta)*x-sin(theta)*y;seed[2*v+1]=sin(theta)*x+cos(theta)*y;
            center[0]+=seed[2*v]/8;center[1]+=seed[2*v+1]/8;
        }
        for (size_t f=0;f<4;f++) { SlimRefineDirection d={{1,0,0},{0,1,0},1};dirs[f]=d; }
        SlimRefineOpts opts;SlimRefine_defaults(&opts);opts.iterations=0;opts.threads=1;opts.strip_pack=0;
        opts.directions=dirs;opts.direction_weight=8;opts.direction_rigid_gauge=1;
        Arena_T arena=Arena_new();SlimRefineStats stats={0};
        int bad=SlimRefine_run_double(arena,xyz,8,faces,4,seed,&opts,result,&stats);
        if (!bad) {
            double co=cos(stats.direction_gauge_radians),si=sin(stats.direction_gauge_radians);
            for (size_t v=0;v<8;v++) {
                double x=seed[2*v]-center[0],y=seed[2*v+1]-center[1];
                if (fabs(result[2*v]-(center[0]+co*x-si*y))>1e-11 ||
                    fabs((v<4?1:-1)*result[2*v+1]-(center[1]+si*x+co*y))>1e-11) bad=1;
            }
        }
        if (bad || !stats.direction_gauge_accepted || fabs(stats.direction_gauge_radians+theta)>1e-10 ||
            fabs(stats.energy_before-stats.energy_after)>1e-10 || stats.direction_angle_after>1e-20 ||
            !(stats.direction_angle_before>.1) || !(stats.objective_after<stats.objective_before)) {
            fprintf(stderr,"[slim direction selftest] mixed-handed common rigid gauge failed\n");failures++;
        }
        Arena_dispose(&arena);
        uint8_t pins[8]={1};opts.fixed_vertices=pins;arena=Arena_new();
        bad=SlimRefine_run_double(arena,xyz,8,faces,4,seed,&opts,result,&stats);
        if (bad || stats.direction_gauge_accepted || result[0]!=seed[0] || result[1]!=seed[1]) {
            fprintf(stderr,"[slim direction selftest] rigid gauge moved an explicit pin\n");failures++;
        }
        Arena_dispose(&arena);
        /* Several distinct quarter-turn assignments create multiple sign
         * changes; compare the exact interval optimum with a full scan. */
        const double phases[4]={-.62,.18,.9,1.4};opts.fixed_vertices=NULL;
        for (size_t f=0;f<4;f++) {
            dirs[f].h[0]=cos(phases[f]);dirs[f].h[1]=sin(phases[f]);
            dirs[f].v[0]=-sin(phases[f]);dirs[f].v[1]=cos(phases[f]);
        }
        arena=Arena_new();bad=SlimRefine_run_double(arena,xyz,8,faces,4,seed,&opts,result,&stats);
        SRSystem system={0};bad|=sr_build_system(arena,xyz,8,faces,4,&opts,1,&system);
        for (int step=-180;!bad && step<=180;step++) {
            double a=(double)step*3.14159265358979323846/720,co=cos(a),si=sin(a),trial[16];
            for (size_t v=0;v<8;v++) {
                double x=seed[2*v]-center[0],y=seed[2*v+1]-center[1];
                trial[2*v]=center[0]+co*x-si*y;
                trial[2*v+1]=(v<4?1:-1)*(center[1]+si*x+co*y);
            }
            SRMeasure measured={0};bad|=sr_measure(&system,trial,&measured);
            if (measured.energy<stats.objective_after-1e-10) bad=1;
        }
        if (bad || !stats.direction_gauge_accepted) {
            fprintf(stderr,"[slim direction selftest] exact common gauge missed a quarter-turn interval\n");failures++;
        }
        Arena_dispose(&arena);
    }
    {
        /* Tangent projection, physical normalization and exact derivatives.
         * The normal components and input arm lengths must have no effect. */
        const float xyz[9]={0,0,0,2,0,0,1,.8f,.6f};
        const int32_t faces[3]={0,1,2};
        double uv[6]={10,-7,12.6,-6.5,10.7,-5.8};
        SlimRefineDirection directions[1]={{{7,-1.8,2.4},{0,4,-2},.75}};
        SlimRefineOpts opts;SlimRefine_defaults(&opts);
        opts.directions=directions;opts.direction_weight=3;
        Arena_T arena=Arena_new();SRSystem system={0};
        SRMeasure before={0},component={0};double product[6];
        int bad=sr_build_system(arena,xyz,3,faces,1,&opts,1,&system);
        if (!bad) bad=sr_assemble(&system,uv,&opts) || sr_measure(&system,uv,&before) ||
                       sr_measure_component(&system,uv,0,&component);
        if (!bad) {
            for (size_t v=0;v<3;v++) sr_block_row(system.diagonal,system.offset,
                system.neighbor,system.adj_edge,system.edge_block,uv,v,product);
            for (size_t k=0;k<6;k++) {
                double saved=uv[k],eps=1e-6;SRMeasure plus={0},minus={0};
                uv[k]=saved+eps;bad|=sr_measure(&system,uv,&plus);
                uv[k]=saved-eps;bad|=sr_measure(&system,uv,&minus);uv[k]=saved;
                double numerical=(plus.energy-minus.energy)/(2*eps);
                double analytical=2*(product[k]-system.rhs[k]);
                if (fabs(numerical-analytical)>2e-7*fmax(1,fabs(analytical))) bad=1;
            }
            if (fabs(before.energy-component.energy)>1e-12 ||
                fabs(before.energy-before.metric_energy-before.direction_energy)>1e-12 ||
                system.direction_faces!=1 || system.direction_area_fraction!=1 ||
                !(before.direction_angle>0)) bad=1;
        }
        if (bad) {
            fprintf(stderr,"[slim direction selftest] physical projection or objective gradient failed\n");failures++;
        }
        Arena_dispose(&arena);
    }
    {
        /* Relabelled and independently reversed arms define the same cross.
         * Exercise both a smooth assignment and its exact 45-degree tie. */
        const float xyz[9]={0,0,0,1,0,0,0,1,0};
        const int32_t faces[3]={0,1,2};
        const double samples[2][6]={{0,0,1.2,.3,-.1,.8},{0,0,1,1,-1,1}};
        for (size_t sample=0;sample<2;sample++) {
            double reference[24]={0};double reference_energy=0,reference_angle=0;
            for (int variant=0;variant<8;variant++) {
                SlimRefineDirection d={0};
                double sh=(variant&2)?-3:2,sv=(variant&4)?-5:4;
                d.h[(variant&1)?1:0]=sh;d.v[(variant&1)?0:1]=sv;d.weight=.6;
                SlimRefineOpts opts;SlimRefine_defaults(&opts);opts.directions=&d;opts.direction_weight=2;
                Arena_T arena=Arena_new();SRSystem system={0};SRMeasure measure={0};double values[24];
                int bad=sr_build_system(arena,xyz,3,faces,1,&opts,1,&system);
                if (!bad) bad=sr_assemble(&system,samples[sample],&opts) || sr_measure(&system,samples[sample],&measure);
                if (!bad) {
                    memcpy(values,system.diagonal,9*sizeof(double));
                    memcpy(values+9,system.edge_block,9*sizeof(double));
                    memcpy(values+18,system.rhs,6*sizeof(double));
                    if (!variant) { memcpy(reference,values,sizeof values);reference_energy=measure.energy;reference_angle=measure.direction_angle; }
                    else {
                        for (size_t k=0;k<24;k++) if (fabs(values[k]-reference[k])>1e-12) bad=1;
                        if (fabs(measure.energy-reference_energy)>1e-12 || fabs(measure.direction_angle-reference_angle)>1e-12) bad=1;
                    }
                }
                if (bad) {
                    fprintf(stderr,"[slim direction selftest] cross symmetry failed sample=%zu variant=%d\n",sample,variant);failures++;
                }
                Arena_dispose(&arena);
            }
        }
    }
    {
        /* Exact isometry is already at the SD lower bound: a direction solve
         * may rotate its free rim, but cannot spend hidden metric energy. */
        const float xyz[12]={0,0,0,1,0,0,1,1,0,0,1,0};
        const int32_t faces[6]={0,1,2,0,2,3};
        SlimRefineDirection directions[2]={{{1,0,0},{0,1,0},1},{{1,0,0},{0,1,0},1}};
        double seed[8],output[8]={0};double angle=.78539816339744830962;
        for (size_t v=0;v<4;v++) {
            seed[2*v]=cos(angle)*xyz[3*v]-sin(angle)*xyz[3*v+1];
            seed[2*v+1]=sin(angle)*xyz[3*v]+cos(angle)*xyz[3*v+1];
        }
        SlimRefineOpts opts;SlimRefine_defaults(&opts);
        opts.directions=directions;opts.direction_weight=8;opts.iterations=150;
        opts.pcg_iterations=100;opts.pcg_tolerance=1e-12;opts.energy_tolerance=1e-12;
        opts.threads=1;opts.strip_pack=0;
        Arena_T arena=Arena_new();SlimRefineStats stats={0};
        int bad=SlimRefine_run_double(arena,xyz,4,faces,2,seed,&opts,output,&stats);
        if (bad || stats.reverted || stats.flips_after || stats.boundary_intersections_after ||
            stats.energy_after>stats.energy_before*(1+1e-10) || !(sqrt(stats.direction_angle_after)<1e-4)) {
            fprintf(stderr,"[slim direction selftest] isometric free-rim turn failed (status=%d reverted=%d SD=%.12g->%.12g angle=%.9g)\n",
                bad,stats.reverted,stats.energy_before,stats.energy_after,sqrt(stats.direction_angle_after));failures++;
        }
        Arena_dispose(&arena);
        /* Unmeasured arms may be absent, but their area cannot be reported as
         * supported. A measured vector normal to the face is an error. */
        directions[1].weight=0;directions[1].h[0]=NAN;directions[1].v[0]=NAN;
        arena=Arena_new();SRSystem system={0};
        bad=sr_build_system(arena,xyz,4,faces,2,&opts,1,&system);
        if (bad || system.direction_faces!=1 || fabs(system.direction_area_fraction-.5)>1e-12) {
            fprintf(stderr,"[slim direction selftest] partial evidence coverage failed\n");failures++;
        }
        Arena_dispose(&arena);
        directions[0].h[0]=0;directions[0].h[2]=1;
        arena=Arena_new();
        if (!sr_build_system(arena,xyz,4,faces,2,&opts,1,&system)) {
            fprintf(stderr,"[slim direction selftest] normal-only material arm accepted\n");failures++;
        }
        Arena_dispose(&arena);
    }
    {
        /* A single connected sheet with a varying material family. Its ideal
         * map has a curved free rim (u=x,v=y-.02*x*x); no rigid turn can align
         * both sides. The initial map is stretched and rotated. */
        enum { NX=9,NY=5,NV=NX*NY,NF=(NX-1)*(NY-1)*2 };
        float xyz[3*NV];int32_t faces[3*NF];double seed[2*NV],output[2*NV]={0},disabled[2*NV]={0};
        SlimRefineDirection directions[NF];size_t f=0;
        for (size_t y=0;y<NY;y++) for (size_t x=0;x<NX;x++) {
            size_t v=y*NX+x;double px=(double)x-4,py=(double)y-2;
            xyz[3*v]=(float)px;xyz[3*v+1]=(float)py;xyz[3*v+2]=0;
            seed[2*v]=4+cos(.3)*1.4*px-sin(.3)*.7*py;
            seed[2*v+1]=-3+sin(.3)*1.4*px+cos(.3)*.7*py;
        }
        for (size_t y=0;y+1<NY;y++) for (size_t x=0;x+1<NX;x++) {
            int32_t a=(int32_t)(y*NX+x),b=a+1,c=a+NX,d=c+1;
            faces[3*f]=a;faces[3*f+1]=b;faces[3*f+2]=d;f++;
            faces[3*f]=a;faces[3*f+1]=d;faces[3*f+2]=c;f++;
        }
        for (f=0;f<NF;f++) {
            /* Each triangle's discrete derivative of x^2 has the cell's
             * midpoint slope, so the target field is exactly integrable. */
            double mid=0;for (size_t k=0;k<3;k++) mid+=xyz[3*(size_t)faces[3*f+k]];
            mid=floor(mid/3)+.5;
            SlimRefineDirection d={{1,.04*mid,0},{0,1,0},1};directions[f]=d;
        }
        SlimRefineOpts opts;SlimRefine_defaults(&opts);
        opts.iterations=150;opts.pcg_iterations=300;opts.pcg_tolerance=1e-11;
        opts.energy_tolerance=1e-11;opts.threads=1;opts.strip_pack=0;
        opts.directions=directions;opts.direction_weight=32;
        Arena_T arena=Arena_new();SlimRefineStats stats={0};
        int bad=SlimRefine_run_double(arena,xyz,NV,faces,NF,seed,&opts,output,&stats);
        double rim_motion=0;
        for (size_t v=0;v<NV;v++) if (v<NX || v>=NV-NX || v%NX==0 || v%NX==NX-1)
            rim_motion+=hypot(output[2*v]-seed[2*v],output[2*v+1]-seed[2*v+1]);
        if (!bad && (stats.reverted || stats.direction_faces!=NF || stats.direction_area_fraction!=1 ||
            stats.flips_after || stats.boundary_intersections_after || !(rim_motion>1) ||
            !(stats.energy_after<stats.energy_before) || !(stats.objective_after<stats.objective_before) ||
            !(sqrt(stats.direction_angle_after)<.02))) bad=1;
        fprintf(stderr,"  SLIM varying direction: SD %.9g -> %.9g, objective %.9g -> %.9g, angular RMS %.6g -> %.6g, rim motion %.6g\n",
            stats.energy_before,stats.energy_after,stats.objective_before,stats.objective_after,
            sqrt(stats.direction_angle_before),sqrt(stats.direction_angle_after),rim_motion);
        Arena_dispose(&arena);
        /* Supplying an evidence array while strength is zero is bitwise
         * equivalent to the legacy solver; no hidden direction enforcement. */
        SlimRefineStats base_stats={0},off_stats={0};opts.direction_weight=0;
        arena=Arena_new();bad|=SlimRefine_run_double(arena,xyz,NV,faces,NF,seed,&opts,disabled,&off_stats);
        Arena_dispose(&arena);opts.directions=NULL;arena=Arena_new();
        bad|=SlimRefine_run_double(arena,xyz,NV,faces,NF,seed,&opts,output,&base_stats);
        if (memcmp(output,disabled,sizeof output) || off_stats.direction_faces || off_stats.direction_energy_after ||
            off_stats.objective_after!=off_stats.energy_after || base_stats.direction_area_fraction) bad=1;
        Arena_dispose(&arena);
        if (bad) { fprintf(stderr,"[slim direction selftest] free varying-field solve or disabled compatibility failed\n");failures++; }
    }
    fprintf(stderr,"  SLIM physical direction objective: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}
