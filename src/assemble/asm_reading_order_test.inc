int AsmReadingOrder_selftest(void)
{
    int fail=0;Arena_T arena=Arena_new();
    AroChart graph[3]={{0}};for(size_t c=0;c<3;c++)graph[c].valid=1;
    AroEdge edges[]={{0,1,1},{1,2,0},{0,2,1}};size_t conflicts=0;
    if(aro_lift(arena,graph,3,edges,3,&conflicts)!=1 || conflicts || graph[2].lift!=1)fail++;
    edges[2].delta=0;
    if(aro_lift(arena,graph,3,edges,3,&conflicts)!=1 || conflicts!=1 || !graph[1].phase_conflict)fail++;
    double origin[3]={0,0,0},direction[3]={1,0,0};AsmAxis *axis=AsmAxis_line(arena,origin,direction,0,256);
    AsmRun run={0};run.arena=arena;run.n_charts=3;run.layer_d_global=10;
    run.charts=ARENA_CALLOC(arena,3,sizeof *run.charts);run.chart_layer_d=ARENA_ALLOC(arena,3*sizeof(double));
    for(size_t c=0;c<3;c++){
        AsmChart *ch=run.charts+c;double r=c==2?130:100+10*(double)c;
        ch->id=(int32_t)c;ch->nv=81;ch->nf=128;ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;run.chart_layer_d[c]=10;
        ch->xyz=ARENA_ALLOC(arena,243*sizeof(float));ch->uv=ARENA_ALLOC(arena,162*sizeof(float));ch->faces=ARENA_ALLOC(arena,384*sizeof(int32_t));
        for(int z=0;z<9;z++)for(int t=0;t<9;t++){
            int v=9*z+t;double theta=-.8+.2*t;
            ch->xyz[3*v]=(float)(32*z);ch->xyz[3*v+1]=(float)(r*cos(theta));ch->xyz[3*v+2]=(float)(r*sin(theta));
            ch->uv[2*v]=(float)(r*theta);ch->uv[2*v+1]=(float)(32*z);
            for(int d=0;d<3;d++)ch->centroid[d]+=ch->xyz[3*v+d]/81.0;
            if(z<8 && t<8){int32_t tris[]={v,v+1,v+10,v,v+10,v+9};memcpy(ch->faces+6*(8*z+t),tris,sizeof tris);}
        }
        for(size_t f=0;f<ch->nf;f++)ch->area3d+=aro_area(ch,f);
    }
    AsmChart saved[3];memcpy(saved,run.charts,sizeof saved);
    float xyz[3][243],uv[3][162];
    for(size_t c=0;c<3;c++){memcpy(xyz[c],run.charts[c].xyz,sizeof xyz[c]);memcpy(uv[c],run.charts[c].uv,sizeof uv[c]);}
    AroModel m={0};if(aro_model(arena,&run,axis,&m) || !m.radial_solved || m.relative_charts!=3 || m.sync.satisfied_relations!=m.sync.eligible_relations)fail++;
    if(m.chart && (fabs((m.chart[1].qlo-m.chart[0].qlo)-1)>1e-5 || fabs((m.chart[2].qlo-m.chart[0].qlo)-3)>1e-5))fail++;
    if(memcmp(saved,run.charts,sizeof saved))fail++;
    for(size_t c=0;c<3;c++)if(memcmp(xyz[c],run.charts[c].xyz,sizeof xyz[c]) || memcmp(uv[c],run.charts[c].uv,sizeof uv[c]))fail++;
    /* No radial comparison is an ordinary unresolved ordering case. */
    AsmRun single=run;single.n_charts=1;
    AroModel one={0};if(aro_model(arena,&single,axis,&one) || one.relative_charts || one.within_charts!=1)fail++;
    AsmRun empty={0};AroModel none={0};if(aro_model(arena,&empty,axis,&none) || none.samples)fail++;
    /* Taper must not reverse the measured handedness of a long helix. */
    { AroChart helix={0};helix.valid=1;helix.count=36;float angles[108];double rr[36],zz[36];size_t used=0;
      for(size_t f=0;f<36;f++){double q=(double)(f%9)*.3,z=(double)(f/9)*40;
        for(int k=0;k<3;k++)angles[3*f+k]=(float)(ARO_TAU*q);rr[f]=100-10*q+.1*z;zz[f]=z;}
      if(aro_sense(arena,&helix,1,1,angles,rr,zz,10,&used)!=-1 || used!=1)fail++; }
    ves_mkdir("output");if(AsmReadingOrder_write(&run,axis,"output/reading_order_selftest"))fail++;
    AsmRibbonKey key[3];
    if(AsmReadingOrder_ribbon_keys(&run,axis,key))fail++;
    else{
        for(size_t c=0;c<3;c++)if(!key[c].source_supported || !(key[c].angular_gradient[0]>0) ||
            fabs(key[c].angular_gradient[1])>1e-8 || !(key[c].turn_hi>key[c].turn_lo))fail++;
        if(fabs(key[1].turn_lo-key[0].turn_lo-1)>.02 || fabs(key[2].turn_lo-key[0].turn_lo-3)>.02)fail++;
    }
    if(memcmp(saved,run.charts,sizeof saved))fail++;
    for(size_t c=0;c<3;c++)if(memcmp(xyz[c],run.charts[c].xyz,sizeof xyz[c]) || memcmp(uv[c],run.charts[c].uv,sizeof uv[c]))fail++;
    run.layer_d_global=0;if(!AsmReadingOrder_ribbon_keys(&run,axis,key))fail++;run.layer_d_global=10;
    fprintf(stderr,"  physical reading order: winding gaps, local pitch, graph cycles, isolated/empty inputs and immutable XYZ/UV: %s (%d failures; %zu/%zu radial relations)\n",
            fail?"FAIL":"ok",fail,m.sync.satisfied_relations,m.sync.eligible_relations);
    Arena_dispose(&arena);return fail;
}
