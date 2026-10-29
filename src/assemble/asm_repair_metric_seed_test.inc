static int amb_selftest(void)
{
    int failures = 0;
    double jacobians[][4] = {{1,0,0,1},{2.2,.4,-.3,.2},{-.3,.4,.2,1.4},{0,-1.5,.4,0},{1,0,0,-1}};
    for (size_t n = 0; n < sizeof jacobians/sizeof jacobians[0]; n++) {
        double residual[4], value = amb_project(jacobians[n], residual), p[4];
        for (int k = 0; k < 4; k++) p[k] = jacobians[n][k]-residual[k];
        double hi = .5*(hypot(p[0]+p[3],p[2]-p[1])+hypot(p[0]-p[3],p[2]+p[1]));
        double det = p[0]*p[3]-p[1]*p[2], lo = det/hi;
        int okay = isfinite(value) && det > 0 && lo >= .9002-1e-12 && hi <= 1.0998+1e-12;
        /* Equal singular values under a reflection have nonunique projections. */
        if (n != 4) for (int k = 0; k < 4; k++) {
            double a[4], b[4], scratch[4]; memcpy(a,jacobians[n],sizeof a); memcpy(b,a,sizeof b);
            a[k] += 1e-6; b[k] -= 1e-6;
            double numeric = (amb_project(a,scratch)-amb_project(b,scratch))/2e-6;
            okay = okay && fabs(numeric-2*residual[k]) < 1e-6*fmax(1,fabs(numeric));
        }
        fprintf(stderr,"  native signed metric projection %zu: %s\n",n,okay?"ok":"FAIL"); failures += !okay;
    }
    for (int variant = 0; variant < 3; variant++) for (int test = 0; test < 5; test++) {
        int n = test == 4 ? 93 : 9; size_t nv = (size_t)n*n, nf = 2*(size_t)(n-1)*(n-1);
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float *xyz = ARENA_ALLOC(arena,3*nv*sizeof(float)), *uv = ARENA_ALLOC(arena,2*nv*sizeof(float));
        int32_t *faces = ARENA_ALLOC(arena,3*nf*sizeof(int32_t));
        for (int y = 0; y < n; y++) for (int x = 0; x < n; x++) {
            int v = y*n+x; xyz[3*v] = uv[2*v] = (float)x; xyz[3*v+1] = uv[2*v+1] = (float)y; xyz[3*v+2] = 0;
            if (x+1 < n && y+1 < n) { int32_t tri[] = {v,v+1,v+n+1,v,v+n+1,v+n}; memcpy(faces+6*(y*(n-1)+x),tri,sizeof tri); }
        }
        AsmChart chart = {0}; chart.id = 0; chart.parent = -1; chart.nv = nv; chart.nf = nf;
        chart.xyz = xyz; chart.uv = uv; chart.faces = faces; chart.placed = 1; chart.placement_state = ASM_PLACE_ROOT;
        run.charts = &chart; run.n_charts = 1;
        AsmField *f = AsmField_new_local(arena,&run,NULL,NULL);
        if (!f) { failures++; Arena_dispose(&arena); continue; }
        size_t v = (size_t)(n/2)*n+n/2; f->uv[2*v] += test == 1 ? 1 : 1.8; f->uv[2*v+1] += .7;
        if (test == 2) for (size_t k = 0; k < 2*nv; k++) f->uv[k] += 1e9;
        if (test == 4) for (size_t k = 0; k < nv; k++) f->uv[2*k] *= .5;
        double *saved = ARENA_ALLOC(arena,2*nv*sizeof(double)); memcpy(saved,f->uv,2*nv*sizeof(double));
        AmrMetric before = amr_metric(f,0,f->uv); AmrResult result = {0};
        int kept = amb_refine(f,0,ves_clock_sec()+(test == 3 ? -1 : 5),variant,&result);
        AmrMetric after = amr_metric(f,0,f->uv);
        int okay = before.bad && (test != 3 ? kept && amr_qualified(&after) :
            !kept && !memcmp(saved,f->uv,2*nv*sizeof(double)) &&
            !strcmp(result.reason,"local_metric_deadline"));
        if(test == 4) okay = okay && result.free_vertices > 8192;
        fprintf(stderr,"  native local metric correction %d/%d: %zu -> %zu bad, %zu steps, %s\n",variant,test,before.bad,after.bad,result.refine_steps,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    return failures;
}
