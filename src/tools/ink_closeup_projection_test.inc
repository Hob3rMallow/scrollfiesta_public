static int ic_projection_check(int ok, const char *name)
{
    fprintf(stderr,"[closeup projection] %s %s\n",ok?"PASS":"FAIL",name);
    return !ok;
}

typedef struct IcProjectionTest {
    IcProjectionLinear cache;
    int mode, ends, outside;
} IcProjectionTest;

static int ic_projection_test_begin(void *context)
{
    IcProjectionTest *t = context;
    if (t->mode == 3) return -1;
    if (t->mode == 4) RAISE(Arena_Failed);
    return ic_projection_begin(&t->cache);
}

static void ic_projection_test_end(void *context)
{
    IcProjectionTest *t = context;
    t->ends++;
    ic_projection_end(&t->cache);
}

static int ic_projection_test_solve(Arena_T arena, void *context, size_t nv, size_t ne,
                                     const double *diagonal, const int32_t *edge_a,
                                     const int32_t *edge_b, const double *blocks,
                                     const double *rhs, double *result)
{
    IcProjectionTest *t = context;
    if (!t->cache.active) t->outside++;
    else {
        /* Force reuse of borrowed scratch between the eight calls. */
        memset(ARENA_ALLOC(arena,100000),0xa5,100000);
        if (t->cache.solves == 1) {
            if (t->mode == 1) return -1;
            if (t->mode == 2) RAISE(Arena_Failed);
        }
    }
    return ic_projection_linear(arena,&t->cache,nv,ne,diagonal,edge_a,edge_b,blocks,rhs,result);
}

static int ic_projection_matrix_selftest(Arena_T arena)
{
    enum { NV = 64, NE = NV-1 };
    double diagonal[3*NV] = {0}, blocks[3*NE] = {0};
    int32_t edge_a[NE], edge_b[NE];
    double rhs[2*NV], out[2*NV] = {0}, reference[2*NV] = {0};
    for (int v = 0; v < NV; v++) {
        diagonal[3*v] = diagonal[3*v+2] = 3;
        rhs[2*v] = 2*v+1; rhs[2*v+1] = 2*v+2;
    }
    for (int e = 0; e < NE; e++) {
        edge_a[e] = e; edge_b[e] = e+1;
        blocks[3*e] = blocks[3*e+2] = -.2;
    }
    IcProjectionLinear c = {0};
    int failures = 0;
    failures += ic_projection_check(!ic_projection_begin(&c) &&
        !ic_projection_linear(arena,&c,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,out) &&
        !ic_slim_linear(arena,NULL,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,reference) &&
        !memcmp(out,reference,sizeof out),"first factor matches one-shot solution exactly");
    SparseFactor_T factor = c.factor;
    Arena_T storage = c.storage;
    failures += ic_projection_check(ic_projection_begin(&c) && c.factor == factor &&
        c.storage == storage && c.active && c.phases == 1,"busy scope preserves its existing owner");
    diagonal[0] += .001;
    failures += ic_projection_check(ic_projection_linear(arena,&c,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,out),
        "changed diagonal refuses stale factor");
    diagonal[0] = 3; blocks[0] = -.21;
    failures += ic_projection_check(ic_projection_linear(arena,&c,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,out),
        "changed block refuses stale factor");
    blocks[0] = -.2; edge_b[0] = 2;
    failures += ic_projection_check(ic_projection_linear(arena,&c,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,out),
        "changed edge refuses stale factor");
    edge_b[0] = 1;
    failures += ic_projection_check(ic_projection_linear(arena,&c,NV-1,NE,diagonal,edge_a,edge_b,blocks,rhs,out) &&
        ic_projection_linear(arena,&c,NV,NE-1,diagonal,edge_a,edge_b,blocks,rhs,out),
        "changed dimensions refuse stale factor");
    rhs[0] = 17;
    failures += ic_projection_check(!ic_projection_linear(arena,&c,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,out) &&
        !ic_slim_linear(arena,NULL,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,reference) &&
        !memcmp(out,reference,sizeof out) && c.factorizations == 1 && c.solves == 2,
        "new RHS reuses the same factor with exact output");
    ic_projection_end(&c);
    failures += ic_projection_check(!c.factor && !c.storage && !c.active &&
        !ic_projection_linear(arena,&c,NV,NE,diagonal,edge_a,edge_b,blocks,rhs,out) &&
        c.factorizations == 1 && c.solves == 2 && !memcmp(out,reference,sizeof out),
        "outside projection the solver remains one-shot with no retained factor");
    failures += ic_projection_check(!ic_projection_begin(&c) &&
        !ic_projection_linear(arena,&c,1,0,diagonal,NULL,NULL,NULL,rhs,out) &&
        fabs(out[0]-rhs[0]/3) < 1e-12 && fabs(out[1]-rhs[1]/3) < 1e-12 && c.factorizations == 2,
        "new scope and zero-edge matrix preserve the small-system backend");
    ic_projection_end(&c);
    return failures;
}

static int ic_projection_selftest(void)
{
    enum { NX = 13, NY = 10, NV = NX*NY, NF = 2*(NX-1)*(NY-1) };
    float xyz[3*NV]; int32_t faces[3*NF];
    double uv[2*NV], reference[2*NV], result[2*NV], sentinel[2*NV];
    SlimRefineDirection arms[NF];
    CloseupAlignOpts options, cached;
    CloseupAlign_defaults(&options);
    options.slim.linear_solve = ic_slim_linear;
    options.slim.iterations = 20; options.slim.threads = 1; options.slim.verbose = 0;
    Arena_T arena = Arena_new_sized(4096);
    int failures = ic_projection_matrix_selftest(arena);
    for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
        int v = y*NX+x;
        xyz[3*v] = (float)x; xyz[3*v+1] = (float)y; xyz[3*v+2] = 0;
    }
    for (int y = 0, f = 0; y < NY-1; y++) for (int x = 0; x < NX-1; x++) {
        int v = y*NX+x;
        faces[3*f] = v; faces[3*f+1] = v+1; faces[3*f+2] = v+NX+1; f++;
        faces[3*f] = v; faces[3*f+1] = v+NX+1; faces[3*f+2] = v+NX; f++;
    }
    for (int mirror = 0; mirror < 2; mirror++) for (int varied = 0; varied < 2; varied++) {
        for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
            int v = y*NX+x;
            double u = 1.08*x+.12*sin(.5*y), w = .95*y+.11*sin(.4*x), angle = .36;
            uv[2*v] = 137+cos(angle)*u-sin(angle)*w;
            uv[2*v+1] = (mirror?-1:1)*(-62+sin(angle)*u+cos(angle)*w);
        }
        memset(arms,0,sizeof arms);
        for (int f = 0; f < NF; f++) {
            double angle = varied?.025*sin(f*.1):0;
            arms[f].h[0] = cos(angle); arms[f].h[1] = sin(angle);
            arms[f].v[0] = -sin(angle); arms[f].v[1] = cos(angle);
            arms[f].weight = varied && !(f%13)?0:2.5;
        }
        IcProjectionLinear cache = {0};
        cached = options; ic_projection_options(&cached,&cache);
        CloseupAlignStats a = {0}, b = {0};
        int r0 = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&options,reference,&a);
        int r1 = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&cached,result,&b);
        failures += ic_projection_check(r0 >= 0 && r0 == r1 && !memcmp(reference,result,sizeof result) &&
            a.selected_seed == b.selected_seed && a.trials == b.trials && a.target_met == b.target_met &&
            a.scale_nonregression == b.scale_nonregression && a.field_feasible == b.field_feasible &&
            a.after.fiber_rms_degrees == b.after.fiber_rms_degrees && a.after.sd_mean == b.after.sd_mean &&
            a.after.bad10_area == b.after.bad10_area && a.after.bad25_area == b.after.bad25_area,
            "complete SLIM output and quality decisions match for both handednesses and partial evidence");
        failures += ic_projection_check(cache.phases == 1 && cache.factorizations == 1 && cache.solves == 8 &&
            !cache.factor && !cache.storage && !cache.active,"one factor serves eight projection solves only");
    }
    IcProjectionTest *test = ARENA_CALLOC(arena,1,sizeof *test);
    cached = options;
    cached.slim.linear_solve = ic_projection_test_solve; cached.slim.linear_context = test;
    cached.projection_begin = ic_projection_test_begin; cached.projection_end = ic_projection_test_end;
    CloseupAlignStats stats = {0}, saved_stats;
    for (int k = 0; k < 2*NV; k++) sentinel[k] = -12345;
    for (int mode = 0; mode <= 4; mode++) {
        memset(test,0,sizeof *test); test->mode = mode;
        memcpy(result,sentinel,sizeof result);
        memset(&stats,0x5a,sizeof stats); memcpy(&saved_stats,&stats,sizeof stats);
        Arena_Mark mark = Arena_save(arena);
        volatile int rc = -99, raised = 0;
        TRY { rc = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&cached,result,&stats); }
        EXCEPT(Arena_Failed) { raised = 1; }
        END_TRY;
        Arena_Mark after = Arena_save(arena);
        failures += ic_projection_check(mark.avail == after.avail && mark.chunk == after.chunk &&
            !test->cache.active && !test->cache.factor && !test->cache.storage &&
            test->ends == (mode < 3 ? 1 : 0),"projection lifetime and borrowed scratch close on every exit");
        if (mode == 2 || mode == 4)
            failures += ic_projection_check(raised && !memcmp(result,sentinel,sizeof result) &&
                !memcmp(&stats,&saved_stats,sizeof stats),"allocation exception propagates with transactional outputs");
        else
            failures += ic_projection_check(!raised && rc >= 0 && test->outside > 0 &&
                test->cache.factorizations == (mode < 3 ? 1 : 0) &&
                test->cache.solves == (mode == 0 ? 8 : mode == 1 ? 1 : 0),
                "projection failure retains independent old-seed SLIM route and never reuses its matrix");
    }
    memset(test,0,sizeof *test);
    for (int i = 0; i < 2; i++)
        failures += ic_projection_check(CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&cached,result,&stats) >= 0,
            "repeated component invocation succeeds");
    failures += ic_projection_check(test->cache.phases == 2 && test->cache.factorizations == 2 &&
        test->cache.solves == 16 && test->ends == 2,"each component owns a separate factor lifetime");
    memset(arms,0,sizeof arms); memset(test,0,sizeof *test);
    failures += ic_projection_check(CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&cached,result,&stats) ==
        CLOSEUP_ALIGN_NO_EVIDENCE && !memcmp(uv,result,sizeof result) && !test->cache.phases,
        "missing evidence never opens a projection scope");
    cached.projection_end = NULL;
    memcpy(result,sentinel,sizeof result); memcpy(&saved_stats,&stats,sizeof stats);
    failures += ic_projection_check(CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&cached,result,&stats) ==
        CLOSEUP_ALIGN_INVALID && !memcmp(result,sentinel,sizeof result) && !memcmp(&stats,&saved_stats,sizeof stats),
        "unpaired hooks refuse before writing outputs");
    Arena_dispose(&arena);
    fprintf(stderr,"[closeup projection] %s: %d failure(s)\n",failures?"FAIL":"PASS",failures);
    return failures;
}
