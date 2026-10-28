/* One unusable proxy vertex must not suppress coherent contraction of the
 * whole strip. Its displacement has an arbitrarily small uniform flip bound;
 * shrinking neighboring alphas instead creates a traveling ring of flips. */
{
    enum {W=65,H=5,NV=W*H,NF=2*(W-1)*(H-1)};
    Arena_T scratch=Arena_new();SRSystem system={0};SlimRefineOpts options;
    float *xyz=ARENA_ALLOC(scratch,3*NV*sizeof(float));
    int32_t *faces=ARENA_ALLOC(scratch,3*NF*sizeof(int32_t));
    double *seed=ARENA_ALLOC(scratch,2*NV*sizeof(double));
    double *destination=ARENA_ALLOC(scratch,2*NV*sizeof(double));
    double *trial=ARENA_ALLOC(scratch,2*NV*sizeof(double));
    double *alpha=ARENA_ALLOC(scratch,NV*sizeof(double));
    uint8_t *damp=ARENA_ALLOC(scratch,NV);
    SREdgeBox *boxes=ARENA_ALLOC(scratch,2*(W+H)*sizeof(SREdgeBox));
    size_t nf=0;
    for (size_t y=0;y<H;y++) for (size_t x=0;x<W;x++) {
        size_t v=y*W+x;
        xyz[3*v]=(float)x;xyz[3*v+1]=(float)y;xyz[3*v+2]=0;
        seed[2*v]=3*(double)x;seed[2*v+1]=.4*(double)y;
        destination[2*v]=(double)x+20;destination[2*v+1]=(double)y;
        if (x+1<W && y+1<H) {
            int32_t a=(int32_t)v,b=a+1,c=a+W,d=c+1;
            faces[3*nf]=a;faces[3*nf+1]=b;faces[3*nf+2]=d;nf++;
            faces[3*nf]=a;faces[3*nf+1]=d;faces[3*nf+2]=c;nf++;
        }
    }
    destination[2*(2*W+W/2)+1]=1e9;
    SlimRefine_defaults(&options);options.threads=1;options.verbose=1;
    SRMeasure before={0},uniform={0},after={0};
    size_t trials=0,moves=0;double amin=0,amean=0,amax=0;
    int ok=!sr_build_system(scratch,xyz,NV,faces,NF,&options,1,&system) &&
           !sr_measure_component(&system,seed,0,&before);
    if (ok) sr_relax_setup(scratch,&system);
    double uniform_alpha=ok?sr_component_flip_step(&system,seed,destination,0):0;
    for (size_t k=0;k<2*NV;k++) trial[k]=seed[k]+uniform_alpha*(destination[k]-seed[k]);
    if (ok) ok=!sr_measure_component(&system,trial,0,&uniform);
    int shrunk=ok && sr_component_local_step(&system,seed,destination,0,&options,
        before.energy,trial,alpha,damp,boxes,&after,&trials,&amin,&amean,&amax);
    int grown=ok && sr_component_growth_step(&system,seed,destination,0,&options,
        before.energy,trial,alpha,boxes,&after,&trials,&moves,&amin,&amean,&amax);
    fprintf(stderr,"[slim growth fixture] shrink=%d growth=%d uniform-alpha=%.9g growth-alpha=%.9g/%.9g E=%.9g/%.9g/%.9g\n",
            shrunk,grown,uniform_alpha,amean,amax,before.energy,uniform.energy,after.energy);
    ok=ok && !shrunk && grown && !after.flips && moves && amean>1000*uniform_alpha &&
       after.energy<uniform.energy && after.energy<before.energy &&
       !sr_component_boundary_intersections(&system,0,trial,boxes,NULL);
    for (size_t f=0;ok && f<NF;f++) ok=UvGuard_interval(seed,trial,faces+3*f,1);
    if (!ok) {fprintf(stderr,"[slim selftest] feasible alpha growth cascade failed\n");failures++;}
    Arena_dispose(&scratch);
}
/* Another edge may have moved many pixels earlier in the same growth sweep.
 * The index must enclose its entire current-to-destination segment, not the
 * one-pixel boxes used by ordinary rim relaxation. */
{
    const double current[8]={0,0,0,2,20,.5,20,.7};
    const double destination[8]={2,0,0,2,.5,.5,.7,.5};
    double trial[8]={2,0,0,2,.5,.5,.7,.5},old[2]={0,0};
    SRBoundaryEdge rim[2]={{0,1,0},{2,3,0}};
    size_t component_offsets[2]={0,2},component_edges[2]={0,1};
    size_t rim_offsets[5]={0,1,2,3,4},rim_edges[4]={0,0,1,1};
    SREdgeBox boxes[2];double prefix[2];SRSystem system={0};
    system.boundary=rim;system.component_boundary_offset=component_offsets;
    system.component_boundary=component_edges;system.rim_offset=rim_offsets;
    system.rim_edge=rim_edges;system.relax_boxes=boxes;system.relax_prefix=prefix;
    sr_growth_rim_index(&system,0,current,destination);
    int bad=sr_relax_rim_safe(&system,0,trial,0,old);
    trial[0]=-.2;
    bad|=!sr_relax_rim_safe(&system,0,trial,0,old);
    if (bad) {fprintf(stderr,"[slim selftest] feasible-growth swept hull failed\n");failures++;}
}
/* Two actual nonadjacent slit-bank edges from the v26 3-iteration replay.
 * Their full endpoints are clear, but uniform interpolation crosses near
 * scale .418. Rebuilding at the half-scale upper bounds must remain clear
 * while making positive progress, rather than accepting the crossing or
 * treating this valid current state as a pre-existing contact. */
{
    const double currents[2][8]={{1367.9131597293404,1490.5602700325071,
        1369.4407956661441,1492.9420768176244,1369.3468259106453,1493.0034529937993,
        1370.8745722046497,1495.3853396199977},
        {1380.134456670552,1509.6886016693547,1381.7170293430447,1512.1669388713783,
         1381.5346095909069,1511.8719049301492,1383.0691529489623,1514.3537762739236}};
    const double shortened[2][8]={{1368.8136785960442,1491.4851694134004,
        1370.0802756089934,1493.9679937935134,1369.6537487136309,1493.1799939096272,
        1370.8919494204706,1495.3949479675607},
        {1380.1344567362034,1509.6886016824174,1381.7170294116995,1512.1669388925761,
         1381.5346097150086,1511.8719050456937,1383.0691530254039,1514.3537764006389}};
    const double scales[2]={.5,1.1920928955078125e-7};
    for (size_t example=0;example<2;example++) {
    const double *current=currents[example],*shortened_trial=shortened[example];
    double destination[8],trial[8],alpha[4],limit[4]={1,1,1,1};
    for (size_t k=0;k<8;k++) destination[k]=current[k]+(shortened_trial[k]-current[k])/scales[example];
    SRBoundaryEdge rim[2]={{0,1,0},{2,3,0}};
    size_t component_offsets[2]={0,2},component_edges[2]={0,1};
    size_t vertex_offsets[2]={0,4},incident_offsets[5]={0,0,0,0,0};
    int32_t vertices[4]={0,1,2,3};uint8_t boundary[4]={1,1,1,1};
    size_t rim_offsets[5]={0,1,2,3,4},rim_edges[4]={0,0,1,1};
    SREdgeBox boxes[2],audit[2];double prefix[2];SRSystem system={0};
    SlimRefineOpts options;SlimRefine_defaults(&options);
    system.boundary=rim;system.component_boundary_offset=component_offsets;
    system.component_boundary=component_edges;system.rim_offset=rim_offsets;
    system.rim_edge=rim_edges;system.relax_boxes=boxes;system.relax_prefix=prefix;
    system.component_vertex_offset=vertex_offsets;system.component_vertex=vertices;
    system.incident_offset=incident_offsets;system.boundary_vertex=boundary;
    int ok=!sr_component_boundary_intersections(&system,0,current,audit,NULL) &&
        !sr_component_boundary_intersections(&system,0,destination,audit,NULL) &&
        sr_component_boundary_intersections(&system,0,shortened_trial,audit,NULL)==1;
    sr_growth_rim_index(&system,0,current,destination);
    size_t face_blocks=0,rim_blocks=0;int sweeps=0;
    size_t moves=sr_growth_candidate(&system,current,destination,0,&options,limit,scales[example],
        trial,alpha,&sweeps,&face_blocks,&rim_blocks);
    double progress=0;
    for (size_t v=0;v<4;v++) {progress+=alpha[v];ok=ok && alpha[v]>=0 && alpha[v]<=scales[example];}
    ok=ok && moves && progress>0 && rim_blocks &&
       !sr_component_boundary_intersections(&system,0,trial,audit,NULL);
    fprintf(stderr,"[slim growth rim fixture %zu] moves=%zu blocks=%zu alpha-sum=%.9g -> %s\n",
        example,moves,rim_blocks,progress,ok?"ok":"FAIL");
    if (!ok) {fprintf(stderr,"[slim selftest] capped rim regrowth failed\n");failures++;}
    }
}
