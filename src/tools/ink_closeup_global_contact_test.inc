/* Same native wrapper/callback and commit helper, with enough triangles that
 * stale inactive material occupies separate tree leaves. The old path is a
 * deliberate negative control, never a production option. */
static int ic_align_global_sequence(int fixed,AsmContactStats *truth,
    int *callback,size_t *leaves)
{
    enum { SIDE=17,PV=SIDE*SIDE,PF=2*(SIDE-1)*(SIDE-1),NV=2*PV,NF=2*PF };
    double uv[2*NV],proposal[2*PV],assembled[2*NV];int32_t faces[3*NF],first[PV],second[PV];
    for (int piece=0;piece<2;piece++) {
        for (int y=0;y<SIDE;y++) for (int x=0;x<SIDE;x++) {
            int v=piece*PV+y*SIDE+x;uv[2*v]=1000*piece+10*x;uv[2*v+1]=10*y;
        }
        int f=piece*PF;
        for (int y=0;y<SIDE-1;y++) for (int x=0;x<SIDE-1;x++) {
            int v=piece*PV+y*SIDE+x;
            faces[3*f]=v;faces[3*f+1]=v+1;faces[3*f+2]=v+SIDE+1;f++;
            faces[3*f]=v;faces[3*f+1]=v+SIDE+1;faces[3*f+2]=v+SIDE;f++;
        }
    }
    for (int v=0;v<PV;v++) { first[v]=v;second[v]=PV+v; }
    MeshBinData mesh={0};mesh.nv=NV;mesh.nf=NF;mesh.faces=faces;
    Arena_T arena=Arena_new();if (!arena) return 1;
    IcAlignGlobalContact *g=ic_align_global_new(arena,&mesh,uv,NULL);
    if (!g) { Arena_dispose(&arena);return 1; }
    *leaves=AsmContacts_leaves(g->query);int bad=*leaves<2;
    g->vertices=first;g->local_vertices=PV;memset(g->moving,1,PF);
    AsmContacts_active_faces(g->query,g->moving);
    for (int v=0;v<PV;v++) { proposal[2*v]=uv[2*v]+100;proposal[2*v+1]=uv[2*v+1]; }
    bad|=ic_align_global_candidate(g,PV,proposal)!=1; /* eventual retained winner */
    for (int v=0;v<PV;v++) proposal[2*v]=uv[2*v]+1000;
    bad|=ic_align_global_candidate(g,PV,proposal)!=0; /* last trial is rejected */
    for (int v=0;v<PV;v++) uv[2*v]+=100;
    if (fixed) bad|=ic_align_global_commit(g)!=0;
    else memcpy(g->trial,uv,sizeof uv); /* v2's original incomplete transition */
    g->vertices=second;memset(g->moving,0,NF);memset(g->moving+PF,1,PF);
    memcpy(assembled,uv,sizeof uv);
    for (int v=0;v<PV;v++) {
        proposal[2*v]=uv[2*(PV+v)]-900;proposal[2*v+1]=uv[2*(PV+v)+1];
        assembled[2*(PV+v)]=proposal[2*v];assembled[2*(PV+v)+1]=proposal[2*v+1];
    }
    *callback=ic_align_global_candidate(g,PV,proposal);
    bad|=memcmp(g->trial,uv,sizeof uv)!=0 || g->unknown!=0;
    AsmContacts *fresh=AsmContacts_new_source(arena,NV,faces,NF,assembled,AsmContacts_material_tolerance(),1,NULL);
    bad|=!fresh || AsmContacts_measure(fresh,assembled,NULL,NULL,NULL,truth)!=0 || !truth->complete || truth->unmeasured;
    Arena_dispose(&arena);return bad;
}

static int ic_align_global_refit_selftest(void)
{
    AsmContactStats stale_truth = {0}, fixed_truth = {0};
    int stale = -1, fixed = -1; size_t first = 0, second = 0;
    int bad = ic_align_global_sequence(0,&stale_truth,&stale,&first) |
              ic_align_global_sequence(1,&fixed_truth,&fixed,&second);
    bad |= stale != 1 || fixed != 0 || !fixed_truth.pairs || fixed_truth.pairs != stale_truth.pairs;
    bad |= first != second || fabs(fixed_truth.area-stale_truth.area) > 1e-9*fmax(1,stale_truth.area);
    fprintf(stderr,"[global guard refit] leaves=%zu old-admits=%d fixed-refuses=%d fresh-pairs=%zu: %s\n",
            first,stale == 1,fixed == 0,fixed_truth.pairs,bad?"FAIL":"PASS");
    return bad != 0;
}

static int ic_align_global_lifetime_selftest(void)
{
    const int32_t faces[12] = {0,1,2,0,2,3,4,5,6,4,6,7}, first[4] = {0,1,2,3};
    double uv[16] = {0,0,20,0,20,20,0,20,40,0,60,0,60,20,40,20}, proposal[8];
    MeshBinData mesh = {0}; mesh.nv = 8; mesh.nf = 4; mesh.faces = (int32_t *)faces;
    Arena_T scratch = Arena_new_sized(4096), owned = Arena_new_sized(4096);
    IcAlignGlobalContact *g = ic_align_global_new(owned,&mesh,uv,NULL);
    int fails = g == NULL;
    if (g) {
        g->vertices = first; g->local_vertices = 4; g->moving[0] = g->moving[1] = 1;
        AsmContacts_active_faces(g->query,g->moving);
        for (int round = 0; round < 6; round++) {
            Arena_Mark mark = Arena_save(scratch);
            double shift = round%2 ? 10 : 26;
            for (size_t k = 0; k < 8; k++) proposal[k] = uv[k]+(k%2 ? 0 : shift);
            fails += ic_align_global_candidate(g,4,proposal) != round%2;
            Arena_restore(scratch,mark);
            for (int k = 0; k < 32; k++) {
                void *poison = ARENA_ALLOC(scratch,4096); memset(poison,0xa5,4096);
            }
            Arena_restore(scratch,mark);
        }
        fails += g->calls != 6 || g->refused != 3 || g->unknown || memcmp(g->trial,uv,sizeof uv) != 0;
        proposal[0] = NAN;
        fails += ic_align_global_candidate(g,4,proposal) != 0 || g->unknown != 1;
        for (size_t k = 0; k < 8; k++) proposal[k] = uv[k]+(k%2 ? 0 : 10);
        fails += ic_align_global_candidate(g,4,proposal) != 1 || ic_align_global_commit(g) != 0;
    }
    Arena_dispose(&owned); Arena_dispose(&scratch);
    fprintf(stderr,"[global guard lifetime] scratch reuse and unknown-trial recovery: %s\n",fails?"FAIL":"PASS");
    return fails != 0;
}

typedef struct IcGlobalVeto { int calls, mode; } IcGlobalVeto;

static int ic_align_global_veto(void *context, size_t nv, const double *uv)
{
    IcGlobalVeto *g = context; (void)nv; (void)uv;
    g->calls++;
    if (g->mode == 2) RAISE(Arena_Failed);
    return !g->mode && g->calls > 1;
}

static int ic_align_global_search_selftest(void)
{
    const float xyz[12] = {0,0,0,20,0,0,20,20,0,0,20,0};
    const int32_t faces[6] = {0,1,2,0,2,3};
    double seed[8], out[8], ordinary[8];
    const SlimRefineDirection arms[2] = {{{1,0,0},{0,1,0},1},{{1,0,0},{0,1,0},1}};
    for (int v = 0; v < 4; v++) {
        seed[2*v] = 100+cos(.35)*xyz[3*v]-sin(.35)*xyz[3*v+1];
        seed[2*v+1] = -70+sin(.35)*xyz[3*v]+cos(.35)*xyz[3*v+1];
    }
    CloseupAlignOpts opts; CloseupAlign_defaults(&opts);
    opts.field_initialize = 0; opts.slim.iterations = 30; opts.slim.linear_solve = ic_slim_linear;
    Arena_T arena = Arena_new();
    CloseupAlignStats baseline = {0}, accepted = {0}, refused = {0};
    int a = CloseupAlign_propose(arena,xyz,4,faces,2,seed,arms,NULL,&opts,ordinary,&baseline);
    int fails = a < 0 || baseline.trials != 1 || !baseline.trial[0].scale_nonregression ||
                !baseline.trial[0].candidate_guard_pass;
    IcGlobalVeto state = {0}; opts.candidate_context = &state; opts.candidate_guard = ic_align_global_veto;
    int b = CloseupAlign_propose(arena,xyz,4,faces,2,seed,arms,NULL,&opts,out,&accepted);
    int continued = b >= 0 && accepted.trials > 1 && accepted.trial[0].scale_nonregression &&
        !accepted.trial[0].candidate_guard_pass && accepted.trial[1].candidate_guard_pass &&
        accepted.trial[1].direction_weight < accepted.trial[0].direction_weight &&
        accepted.after.fiber_rms_degrees < accepted.before.fiber_rms_degrees;
    fails += !continued;
    state.calls = 0; state.mode = 1;
    int c = CloseupAlign_propose(arena,xyz,4,faces,2,seed,arms,NULL,&opts,out,&refused);
    int retained = c == CLOSEUP_ALIGN_UNCHANGED && refused.trials == 4 && state.calls == 4 &&
                   !memcmp(out,seed,sizeof seed);
    fails += !retained;
    /* The physical wrapper must compose an upstream filter with its material
     * check, and its allocation failure must not commit partial output. */
    MeshBinData mesh = {0}; mesh.nv = 4; mesh.nf = 2;
    mesh.verts = (float *)xyz; mesh.faces = (int32_t *)faces; mesh.uv64 = seed;
    SlimRefineStats stats = {0}, saved_stats; IcPhysicalAlignStats info = {0}, saved_info;
    state.calls = 0;
    int wrapped = ic_reparameterize_physical(&mesh,arms,NULL,&opts,out,&stats,&info,NULL);
    fails += wrapped || !state.calls || memcmp(out,seed,sizeof seed) ||
             info.global_contact_refusals != (size_t)state.calls;
    for (int k = 0; k < 8; k++) out[k] = ordinary[k] = -12345;
    memset(&stats,0x5a,sizeof stats); memcpy(&saved_stats,&stats,sizeof stats);
    memset(&info,0x5a,sizeof info); memcpy(&saved_info,&info,sizeof info);
    state.calls = 0; state.mode = 2;
    wrapped = ic_reparameterize_physical(&mesh,arms,NULL,&opts,out,&stats,&info,NULL);
    fails += !wrapped || !state.calls || memcmp(out,ordinary,sizeof out) ||
             memcmp(&stats,&saved_stats,sizeof stats) || memcmp(&info,&saved_info,sizeof info);
    Arena_dispose(&arena);
    fprintf(stderr,"[global guard search] continue=%d exact-fallback=%d hook-composition and exception transaction: %s\n",
            continued,retained,fails?"FAIL":"PASS");
    return fails != 0;
}

static int ic_align_global_seed_selftest(void)
{
    float xyz[24] = {0,0,0,20,0,0,20,20,0,0,20,0,0,0,10,20,0,10,20,20,10,0,20,10};
    int32_t faces[12] = {0,1,2,0,2,3,4,5,6,4,6,7}, aliases[8] = {0,1,2,3,4,5,6,8};
    double uv[16] = {0,0,20,0,20,20,0,20,14,0,34,0,34,20,14,20}, out[16], sentinel[16];
    MeshBinData mesh = {0}; mesh.nv = 8; mesh.nf = 4; mesh.verts = xyz; mesh.faces = faces; mesh.uv64 = uv;
    CloseupAlignOpts opts; CloseupAlign_defaults(&opts); opts.slim.linear_solve = ic_slim_linear;
    for (int k = 0; k < 16; k++) out[k] = sentinel[k] = -12345;
    SlimRefineStats stats, saved_stats; IcPhysicalAlignStats info, saved_info;
    memset(&stats,0x5a,sizeof stats); memcpy(&saved_stats,&stats,sizeof stats);
    memset(&info,0x5a,sizeof info); memcpy(&saved_info,&info,sizeof info);
    int rc = ic_reparameterize_physical(&mesh,NULL,NULL,&opts,out,&stats,&info,NULL);
    int fails = !rc || memcmp(out,sentinel,sizeof out) || memcmp(&stats,&saved_stats,sizeof stats) ||
                 memcmp(&info,&saved_info,sizeof info);
    for (int v = 4; v < 8; v++) uv[2*v] += 30;
    rc = ic_reparameterize_physical_material(&mesh,NULL,NULL,&opts,out,&stats,&info,NULL,aliases,1);
    fails += !rc || memcmp(out,sentinel,sizeof out) || memcmp(&stats,&saved_stats,sizeof stats);
    aliases[7] = 7;
    rc = ic_reparameterize_physical_material(&mesh,NULL,NULL,&opts,out,&stats,&info,NULL,aliases,1);
    fails += rc || memcmp(out,uv,sizeof out) || info.unknown != 2 || info.global_contact_calls;
    fprintf(stderr,"[global guard seed] unsafe seed/invalid aliases refuse atomically; unknown evidence preserved: %s\n",
            fails?"FAIL":"PASS");
    return fails != 0;
}

static int ic_align_global_selftest(void)
{
    int fails = ic_align_global_refit_selftest();
    fails += ic_align_global_lifetime_selftest();
    fails += ic_align_global_search_selftest();
    fails += ic_align_global_seed_selftest();
    fprintf(stderr,"[global component guard] %s: %d failed groups\n",fails?"FAIL":"PASS",fails);
    return fails != 0;
}
