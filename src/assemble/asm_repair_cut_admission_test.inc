static int api_cut_owner_selftest(void)
{
    enum {N=8}; Arena_T arena = Arena_new(), temporary = Arena_new();
    float xyz[3*(N+1)] = {0}, uv[2*(N+1)] = {0}; int32_t vid[N+1], faces[3*N];
    xyz[2] = 1;
    for (int v = 0; v <= N; v++) vid[v] = v;
    for (int k = 0; k < N; k++) {
        double angle = 2*acos(-1.)*k/N;
        xyz[3*(k+1)] = uv[2*(k+1)] = (float)(.25*cos(angle));
        xyz[3*(k+1)+1] = uv[2*(k+1)+1] = (float)(.25*sin(angle));
        faces[3*k] = 0; faces[3*k+1] = k+1; faces[3*k+2] = (k+1)%N+1;
    }
    AsmChart ch = {0}; ch.id = 0; ch.parent = -1; ch.nv = N+1; ch.nf = N;
    ch.xyz = xyz; ch.uv = uv; ch.vid = vid; ch.faces = faces;
    AsmRun original = {0}; original.arena = arena; original.charts = &ch; original.n_charts = 1;
    AsmRun trial = original; trial.arena = temporary;
    trial.charts = ARENA_ALLOC(temporary,sizeof ch); *trial.charts = ch;
    trial.charts->placed = 1; trial.charts->placement_state = ASM_PLACE_SEAM;
    AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1; opts.deadline = ves_clock_sec()+3;
    AprWorkspace work = {0}; opts.workspace = &work; const char *reason = NULL; uint8_t core = 1;
    int okay = apa_intrinsic(&trial,&core,NULL,&opts,&reason,NULL); AsmRun staged = {0};
    okay = okay && trial.charts->nv > ch.nv && trial.n_rels && api_cut_stage(&original,&trial,&core,&staged);
    if (okay) {
        Arena_dispose(&temporary);
        AsmField *field = AsmField_new(arena,&staged); AsmContactStats contact = {0};
        okay = field && apa_chart_metric(field,0) && amb_cuts_pass(field) &&
            !AsmContacts_measure(field->contacts,field->uv,NULL,NULL,NULL,&contact) && !contact.pairs &&
            apr_cut_preserved(arena,&ch,&staged,0) && !ch.placed && original.n_rels == 0 && original.charts == &ch;
        if (field) {
            float saved = staged.charts->xyz[3*ch.nv]; staged.charts->xyz[3*ch.nv] += .01f;
            okay = okay && !apr_cut_preserved(arena,&ch,&staged,0); staged.charts->xyz[3*ch.nv] = saved;
            unsigned continuity = staged.rels[0].continuity; staged.rels[0].continuity = 0;
            okay = okay && !apr_cut_preserved(arena,&ch,&staged,0); staged.rels[0].continuity = continuity;
        }
    }
    fprintf(stderr,"  native automatic cut recovery, transaction ownership and independent preservation: %s (%s)\n",okay?"ok":"FAIL",reason?reason:"qualified");
    if (temporary) Arena_dispose(&temporary); Arena_dispose(&arena); return !okay;
}

static int api_cut_overlay_selftest(void)
{
    int failures = 0;
    for (int batched = 0; batched < 2; batched++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[] = {0,0,0,12,0,0,12,12,0,0,12,0,0,0,0,12,12,0};
        float uv[] = {0,0,12,0,12,12,1,11,1,-1,13,11};
        int32_t faces[] = {0,1,2,4,5,3}, vid[] = {0,1,2,3,0,2};
        float other_xyz[] = {100,100,0,101,100,0,100,101,0}, other_uv[] = {100,100,101,100,100,101}; int32_t other_faces[] = {0,1,2};
        AsmChart charts[2] = {{0}}; run.charts = charts; run.n_charts = 2;
        charts[0].id = 0; charts[0].parent = -1; charts[0].nv = 6; charts[0].nf = 2;
        charts[0].xyz = xyz; charts[0].uv = uv; charts[0].vid = vid; charts[0].faces = faces;
        charts[0].placed = 1; charts[0].placement_state = ASM_PLACE_ROOT;
        charts[1] = charts[0]; charts[1].id = 1; charts[1].nv = 3; charts[1].nf = 1;
        charts[1].xyz = other_xyz; charts[1].uv = other_uv; charts[1].faces = other_faces;
        AsmRelation r = {0}; r.flags = ASM_REL_WEAK; r.corr_count = r.n_corr = 2;
        r.continuity = ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_REQUIRED | ASM_CONT_INTRACHART_CUT;
        AsmCorr corr[2] = {{0}};
        for (int k = 0; k < 2; k++) {
            corr[k].va = k?2:0; corr[k].vb = k?5:4; corr[k].valid = 3;
            corr[k].trim_a.valid = corr[k].trim_b.valid = 1; corr[k].trim_b.face = 1;
        }
        run.rels = &r; run.n_rels = 1; run.corr = corr; run.n_corr = 2;
        ApaOverlay overlay = {0}; overlay.arena = arena; overlay.batched = batched;
        uint8_t first[] = {1,0}, second[] = {0,1};
        int okay = !apa_overlay_add(&overlay,&run,first) && !apa_overlay_add(&overlay,&run,second);
        AsmContacts *tree = batched ? overlay.level[1]->tree : overlay.tree;
        double *encoded = batched ? overlay.level[1]->uv : overlay.uv; AsmContactStats contact = {0};
        okay = okay && tree && !AsmContacts_measure(tree,encoded,NULL,NULL,NULL,&contact) && contact.pairs == 1;
        fprintf(stderr,"  native cut material adjacency survives overlay merge %d: %zu pairs, %s\n",batched,contact.pairs,okay?"ok":"FAIL");
        failures += !okay; apa_overlay_clear(&overlay); Arena_dispose(&arena);
    }
    return failures;
}

static int api_cut_selftest(void)
{ return api_cut_owner_selftest()+api_cut_overlay_selftest(); }
