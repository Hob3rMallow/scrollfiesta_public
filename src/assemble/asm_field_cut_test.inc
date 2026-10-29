static int aff_cut_material_selftest(void)
{
    int failures = 0;
    for (int test = 0; test < 4; test++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[] = {0,0,0, 12,0,0, 12,12,0, 0,12,0, 0,0,0, 12,12,0};
        float uv[] = {0,0, 12,0, 12,12, 1,11, 1,-1, 13,11};
        int32_t faces[] = {0,1,2, 4,5,3}, vid[] = {0,1,2,3,0,2};
        AsmChart ch = {0}; ch.id = 0; ch.parent = -1; ch.nv = 6; ch.nf = 2;
        ch.xyz = xyz; ch.uv = uv; ch.vid = vid; ch.faces = faces; ch.placed = 1; ch.placement_state = ASM_PLACE_ROOT;
        if (test == 1) { ch.flags |= ASM_CHART_MIRROR; ch.pose_x = 1e9; ch.pose_y = -1e9; }
        AsmRelation r = {0}; r.flags = ASM_REL_WEAK; r.corr_count = r.n_corr = 2;
        r.continuity = ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_REQUIRED | ASM_CONT_INTRACHART_CUT;
        AsmCorr corr[2] = {{0}};
        for (int k = 0; k < 2; k++) {
            corr[k].va = k?2:0; corr[k].vb = k?5:4; corr[k].valid = 3;
            corr[k].trim_a.valid = corr[k].trim_b.valid = 1; corr[k].trim_b.face = 1;
        }
        if (test == 2) vid[4] = 3;
        /* A torn cut is no seam, but the fold across it still counts. */
        if (test == 3) r.flags |= ASM_REL_TORN;
        run.charts = &ch; run.n_charts = 1; run.rels = &r; run.n_rels = 1; run.corr = corr; run.n_corr = 2;
        AsmField *f = AsmField_new_local(arena,&run,NULL,NULL); int okay = f != NULL;
        AsmContactStats raw = {0}, original = {0};
        if (f) {
            int status = AsmField_build_contacts(f);
            if (test == 2) okay = status != 0;
            else if (test == 3) {
                AsmContacts *split = AsmContacts_new(arena,f->nv,f->faces,f->nf,f->uv,AsmContacts_material_tolerance());
                okay = status == 0 && split && f->n_seams == 0 &&
                    !AsmContacts_measure(split,f->uv,NULL,NULL,NULL,&raw) && !raw.pairs &&
                    !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&original) && original.pairs == 1 && original.area > 0;
            }
            else {
                AsmContacts *split = AsmContacts_new(arena,f->nv,f->faces,f->nf,f->uv,AsmContacts_material_tolerance());
                AsmSeamMeasure seam; AsmField_seam(f,0,f->uv,&seam);
                okay = status == 0 && split && seam.pass &&
                    !AsmContacts_measure(split,f->uv,NULL,NULL,NULL,&raw) && !raw.pairs &&
                    !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&original) && original.pairs == 1 && original.area > 0;
            }
        }
        fprintf(stderr,"  native UV cut preserves original material boundary %d: %zu split / %zu original pairs, %s\n",
            test,raw.pairs,original.pairs,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    return failures;
}

static int aff_cut_selftest(void)
{
    int failures = aff_cut_material_selftest();
    for (int test = 0; test < 17; test++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[] = {0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,0, 1,1,0};
        float uv[] = {0,0, 1,0, 1,1, 0,1, 0,0, 1,1};
        int32_t faces[] = {0,1,2, 4,5,3}, vid[] = {0,1,2,3,0,2};
        AsmChart ch = {0}; ch.id = 0; ch.parent = -1; ch.nv = 6; ch.nf = 2;
        ch.xyz = xyz; ch.uv = uv; ch.vid = vid; ch.faces = faces; ch.placed = 1; ch.placement_state = ASM_PLACE_ROOT;
        AsmRelation r = {0}; r.flags = ASM_REL_WEAK; r.corr_count = r.n_corr = 2;
        r.continuity = ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_REQUIRED | ASM_CONT_INTRACHART_CUT;
        AsmCorr corr[3] = {{0}};
        for (int k = 0; k < 2; k++) {
            corr[k].va = k?2:0; corr[k].vb = k?5:4; corr[k].valid = 3;
            corr[k].trim_a.valid = corr[k].trim_b.valid = 1; corr[k].trim_b.face = 1;
        }
        if (test == 1) r.continuity &= ~ASM_CONT_INTRACHART_CUT;
        if (test == 2) vid[4] = 3;
        if (test == 3) xyz[12] = .01f;
        if (test == 4) corr[0].gap_a[0] = .01f;
        if (test == 5) corr[0].trim_b.face = 0;
        if (test == 6) { faces[4] = 3; faces[5] = 5; }
        if (test == 7) r.corr_count = 1;
        if (test == 8) r.continuity &= ~ASM_CONT_REQUIRED;
        if (test == 9) uv[8] = 4.01f;
        if (test == 10) uv[8] = 4;
        if (test == 11) r.flags |= ASM_REL_CONTACT;
        if (test == 12) r.flags |= ASM_REL_DROPPED;
        if (test == 13) r.corr_count = 3;
        if (test == 14) corr[0].trim_a.inverse[0] = NAN;
        if (test == 15) { ch.flags |= ASM_CHART_MIRROR; ch.pose_x = 1e9; ch.pose_y = -1e9; }
        if (test == 16) vid[0] = vid[4] = -1;
        run.charts = &ch; run.n_charts = 1; run.rels = &r; run.n_rels = 1; run.corr = corr; run.n_corr = 3;
        AsmField *f = AsmField_new_local(arena,&run,NULL,NULL); AsmSeamMeasure m = {0};
        if (f && f->n_seams == 1) AsmField_seam(f,0,f->uv,&m);
        int expected = test == 0 || test == 10 || test == 15;
        int okay = test == 1 ? f == NULL : f && f->n_seams == 1 && f->seams[0].required && m.pass == expected;
        AsmSeamMeasure direct={0};int direct_rc=AsmField_measure_source_seams(&run,&direct);
        if ((test==1 && direct_rc!=-1) || (test!=1 && (direct_rc || memcmp(&m,&direct,sizeof m)))) okay=0;
        if (expected) okay = okay && m.observations == 2 && fabs(m.mass-sqrt(2.)) < 1e-12;
        fprintf(stderr,"  native original-edge cut %d: %s\n",test,okay?"ok":"FAIL"); failures += !okay;
        Arena_dispose(&arena);
    }
    return failures;
}
