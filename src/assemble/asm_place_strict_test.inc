static int aps_selftest(void)
{
    int failures = 0;
    for (int test = 0; test < 18; test++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 12;
        float xyz[3][243], uv[3][162], original_uv[3][162], normals[243];
        int32_t faces[3][384]; uint8_t boundary[81]; double current[162];
        for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++) {
            int v = j*9 + u; boundary[v] = j == 0 || j == 8 || u == 0 || u == 8;
            normals[3*v] = normals[3*v+1] = 0; normals[3*v+2] = 1;
            for (int c = 0; c < 3; c++) {
                xyz[c][3*v] = (float)(j*16); xyz[c][3*v+1] = (float)(u*5 + c*44); xyz[c][3*v+2] = 0;
                uv[c][2*v] = (float)(u*5); uv[c][2*v+1] = (float)(j*16);
                if (j < 8 && u < 8) {
                    int32_t f[6] = {v,v+1,v+10,v,v+10,v+9}; memcpy(faces[c] + (j*8+u)*6, f, sizeof f);
                }
            }
            current[2*v] = u*5 + 44; current[2*v+1] = j*16;
        }
        if (test == 2) for (int v = 0; v < 81; v++) uv[1][2*v] *= 2;
        for (int c = 0; c < 3; c++) {
            AsmChart ch = {0}; ch.id = c; ch.component = c == 1 ? 1 : 0; ch.placed = 1;
            ch.nv = 81; ch.nf = 128; ch.xyz = xyz[c]; ch.uv = uv[c]; ch.faces = faces[c]; ch.nrm = normals; ch.boundary = boundary;
            ch.area3d = 5120; ch.pose_x = c == 0 ? 0 : c == 1 ? 500 : 200;
            ch.placement_state = c == 1 ? ASM_PLACE_NONE : ASM_PLACE_ROOT;
            if (c == 2 && (test == 3 || test == 11 || test == 14)) ch.placement_state = ASM_PLACE_NONE;
            if (c == 2 && test == 1) ch.pose_x = 44;
            if (c == 2 && test == 8) ch.placed_uv = current;
            AsmRun_push_chart(&run, &ch);
        }
        for (int i = 0; i < (test >= 3 && test <= 6 ? 2 : 1); i++) {
            AsmRelation rel = {0}; rel.a = i; rel.b = i+1; rel.flags = i ? ASM_REL_WEAK : 0;
            rel.rms = .1; rel.corr_first = (int32_t)run.n_corr; rel.corr_count = 9; rel.n_corr = 9;
            for (int j = 0; j < 9; j++) { AsmCorr c = {j*9+8,j*9}; AsmRun_push_corr(&run, &c); }
            AsmRun_push_rel(&run, &rel);
        }
        int good = !AsmContinuity_prepare(&run);
        if (test == 12) memcpy(faces[1] + 3, faces[1], 3*sizeof(int32_t));
        if (test == 3 || test == 4) run.rels[1].continuity |= ASM_CONT_REQUIRED;
        AsmChart proposal[2]; proposal[0] = run.charts[1]; proposal[0].pose_x = 44; proposal[0].placement_state = ASM_PLACE_SEAM;
        if (test == 10) proposal[0].pose_x = NAN;
        int32_t groups[3] = {0,1,0}; uint8_t selected[2] = {1,0};
        if (test == 6) selected[1] = 1;
        if (test == 11) { groups[2] = 1; proposal[1] = run.charts[2]; proposal[1].pose_x = 44; proposal[1].placement_state = ASM_PLACE_LAYER; }
        if (test == 13) { for (int c = 0; c < 3; c++) run.charts[c].pose_x += 1e9; proposal[0].pose_x += 1e9; }
        if (test == 16) {
            for (int c = 0; c < 3; c++) { run.charts[c].flags |= ASM_CHART_MIRROR; run.charts[c].pose_x = -run.charts[c].pose_x; }
            proposal[0].flags |= ASM_CHART_MIRROR; proposal[0].pose_x = -proposal[0].pose_x;
        }
        AplCtx x = {0}; x.run = &run; x.orig_comp = groups; x.proposal = proposal; x.nproposal = test == 11 ? 2 : 1; x.seam_selected = selected;
        x.strict = aps_new(run.arena, &run, NULL);
        if (test == 7) x.strict->field_cap = 1;
        if (test == 9) x.strict->deadline = ves_clock_sec() - 1;
        if (test == 14) { run.charts[2].placement_state = ASM_PLACE_LAYER; run.charts[2].pose_x = 44; good &= aps_box(run.charts + 2, x.strict->boxes + 2); }
        if (test == 15) x.strict->walk_budget = 1;
        if (test == 17) x.strict->candidate_cap = 1;
        AsmChart before[3]; AsmRelation rels[2];
        memcpy(before, run.charts, sizeof before); memcpy(rels, run.rels, run.n_rels*sizeof *rels); memcpy(original_uv, uv, sizeof uv);
        ApsResult local, full; int passed = aps_check(&x, 1, 0, &local);
        int expected = test == 0 || test == 3 || test == 5 || test == 13 || test == 16;
        good &= passed == expected;
        if (test != 7 && test != 9 && test != 15 && test != 17) {
            int oracle = aps_check(&x, 1, 1, &full);
            good &= oracle == passed && local.bad_faces == full.bad_faces && local.failed_seams == full.failed_seams &&
                local.contacts.complete == full.contacts.complete && local.contacts.pairs == full.contacts.pairs &&
                fabs(local.contacts.area - full.contacts.area) < 1e-8;
        }
        if (test == 15) good &= !local.contacts.complete && !strcmp(local.reason, "strict_contact_incomplete");
        good &= !memcmp(before, run.charts, sizeof before) && !memcmp(rels, run.rels, run.n_rels*sizeof *rels) && !memcmp(original_uv, uv, sizeof uv);
        aps_require_selected(&run, selected);
        good &= (run.rels[0].continuity & ASM_CONT_REQUIRED) != 0;
        if (run.n_rels > 1) good &= !!(run.rels[1].continuity & ASM_CONT_REQUIRED) == (test == 3 || test == 4 || test == 6);
        if (!good) { fprintf(stderr, "  strict placement selftest FAIL: case %d, expected %d got %d (%s), bad %zu seams %zu contacts %zu\n",
            test, expected, passed, local.reason, local.bad_faces, local.failed_seams, local.contacts.pairs); failures++; }
        Arena_dispose(&run.arena);
    }
    fprintf(stderr, "  strict placement certificates: %d failures across 18 controls\n", failures);
    return failures;
}
