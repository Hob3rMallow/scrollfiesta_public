/* The same planar material beyond the old free-vertex cutoff must be
 * recoverable. Memory/deadline refusals preserve the caller's cached seed;
 * a provably impossible closed source fan is classified before fitting. */
static int api_selftest(void)
{
    enum {N=97, NV=N*N, NF=2*(N-1)*(N-1)}; int failures = 0;
    for (int test = 0; test < 4; test++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float *xyz = ARENA_CALLOC(arena,3*NV,sizeof(float)), *uv = ARENA_ALLOC(arena,2*NV*sizeof(float));
        int32_t *faces = ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
            int v = y*N+x; xyz[3*v] = (float)x; xyz[3*v+1] = (float)y;
            uv[2*v] = .85f*x; uv[2*v+1] = .85f*y;
            if (x+1 < N && y+1 < N) {
                int32_t tri[] = {v,v+1,v+N+1,v,v+N+1,v+N};
                memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);
            }
        }
        AsmChart chart = {0}; chart.id = 0; chart.parent = -1; chart.nv = NV; chart.nf = NF;
        chart.xyz = xyz; chart.uv = uv; chart.faces = faces; chart.placed = 1; chart.placement_state = ASM_PLACE_ROOT;
        if (test == 1) chart.flags |= ASM_CHART_MIRROR;
        run.charts = &chart; run.n_charts = 1;
        float *saved = ARENA_ALLOC(arena,2*NV*sizeof(float)); memcpy(saved,uv,2*NV*sizeof(float));
        AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1; opts.deadline = ves_clock_sec()+5;
        if (test == 2) opts.patch_bytes = 1;
        if (test == 3) opts.deadline = ves_clock_sec()-1;
        uint8_t core[] = {1}; size_t oversize = 0; const char *reason = "qualified";
        int recovered = apa_intrinsic(&run,core,&oversize,&opts,&reason,NULL);
        int okay = recovered == (test < 2) && !memcmp(saved,uv,2*NV*sizeof(float));
        if (recovered) {
            AsmField *field = AsmField_new_local(arena,&run,NULL,NULL);
            okay = okay && chart.uv != uv && field && apa_chart_metric(field,0) && chart.nv == NV && chart.nf == NF && chart.xyz == xyz && chart.faces == faces;
        } else okay = okay && chart.uv == uv && !strcmp(reason,test == 2 ? "workspace_memory_budget" : "intrinsic_recovery_deadline");
        fprintf(stderr,"  native bounded intrinsic grid %d (%d vertices): %s, %s\n",test,NV,reason,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[] = {0,0,3, -1,-1,0, 1,-1,0, 1,1,0, -1,1,0};
        float uv[] = {0,0, -1,-1, 1,-1, 1,1, -1,1}; int32_t faces[] = {0,1,2,0,2,3,0,3,4,0,4,1};
        AsmChart chart = {0}; chart.id = 0; chart.parent = -1; chart.nv = 5; chart.nf = 4;
        chart.xyz = xyz; chart.uv = uv; chart.faces = faces; chart.placed = 1; chart.placement_state = ASM_PLACE_ROOT;
        run.charts = &chart; run.n_charts = 1;
        AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1; opts.deadline = ves_clock_sec()+1;
        uint8_t core[] = {1}; size_t oversize = 0; const char *reason = "unset";
        int recovered = apa_intrinsic(&run,core,&oversize,&opts,&reason,NULL);
        int okay = !recovered && !strcmp(reason,"intrinsic_cut_source_identity_unavailable") && chart.uv == uv && !chart.placed_uv;
        fprintf(stderr,"  native obstructed intrinsic fan: %s, %s\n",reason,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    {
        /* Excluding a single turn is not sufficient when the admitted
         * material band allows this tiny double-winding fan. All positive
         * windings must be excluded before bypassing the actual fitter. */
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[27] = {0}, uv[18] = {0}; int32_t faces[24];
        for (int k = 0; k < 8; k++) {
            double angle = .5*acos(-1.)*k;
            xyz[3*(k+1)] = (float)(.25*cos(angle)); xyz[3*(k+1)+1] = (float)(.25*sin(angle));
            uv[2*(k+1)] = .5f*xyz[3*(k+1)]; uv[2*(k+1)+1] = .5f*xyz[3*(k+1)+1];
            faces[3*k] = 0; faces[3*k+1] = k+1; faces[3*k+2] = (k+1)%8+1;
        }
        AsmChart chart = {0}; chart.id = 0; chart.parent = -1; chart.nv = 9; chart.nf = 8;
        chart.xyz = xyz; chart.uv = uv; chart.faces = faces; chart.placed = 1; chart.placement_state = ASM_PLACE_ROOT;
        run.charts = &chart; run.n_charts = 1; AmfStats source = amf_chart(&chart,NULL);
        AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1; opts.deadline = ves_clock_sec()+1;
        uint8_t core[] = {1}; size_t oversize = 0; const char *reason = "qualified";
        int recovered = apa_intrinsic(&run,core,&oversize,&opts,&reason,NULL);
        AsmField *field = recovered ? AsmField_new(arena,&run) : NULL; AsmContactStats contacts = {0};
        int okay = source.io && source.obstruction_stars == 1 && !source.orientation_obstructions && recovered &&
            field && apa_chart_metric(field,0) && !AsmContacts_measure(field->contacts,field->uv,NULL,NULL,NULL,&contacts) && !contacts.pairs;
        fprintf(stderr,"  native intrinsic multiple-winding allowance: %s, %s\n",reason,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    return failures;
}
