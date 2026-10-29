/* Source topology is a private proposal until ordinary admission proves
 * placement, every incident obligation and complete material contacts.
 * The opened chart is re-flattened by the metric-band refinement with
 * precondition `first` (1: the chart metric alone, 2: metric, cut gaps and
 * contacts jointly), then by the joint refinement if the first did not
 * qualify; each result must pass the same float-metric, cut-gap and contact
 * certificate. Admission starts at 2 and never tears. With `tear` the cut
 * seams are recorded tears (ASM_REL_TORN, no obligations): the slit may open
 * as a cone's development does, and only the metric and contacts remain.
 * `prepared`, when given, reports whether the cut itself was made. */
static int api_cut_recover_from(AsmRun *trial, size_t chart, const AsmRepairOpts *opts, int first, int tear,
                                const char **reason, FILE *trace, int *prepared)
{
    double start = ves_clock_sec(); AsmRepairOpts budget = *opts;
    if (!(budget.deadline > 0)) budget.deadline = start+budget.patch_seconds;
    if (budget.admission_solver_deadline > 0 &&
        (!(budget.deadline > 0) || budget.deadline > budget.admission_solver_deadline)) budget.deadline = budget.admission_solver_deadline;
    AsmRun opened = {0}; AcutResult cut = {0}; AmrResult result = {0};
    const char *stop = "intrinsic_cut_failed"; int okay = 0; AsmContactStats contacts = {0};
    if (first < 1 || first > 2) goto done;
    if (chart < trial->n_charts && !trial->charts[chart].vid) { stop = "intrinsic_cut_source_identity_unavailable"; goto done; }
    if (prepared) *prepared = 0;
    if (!acut_make_mode_budget(trial,chart,SIZE_MAX,trial->arena,&opened,&cut,0,&budget)) { stop = cut.reason; goto done; }
    if (prepared) *prepared = 1;
    if (tear) for (size_t i = trial->n_rels; i < opened.n_rels; i++) opened.rels[i].flags |= ASM_REL_TORN;
    for (int precondition = first; precondition <= 2 && !okay; precondition++) {
        if (precondition > first && apr_expired(&budget)) break;
        Arena_T scratch = Arena_new(); AsmRun one = {0}; one.arena = scratch; one.n_charts = 1;
        AsmChart local = opened.charts[chart]; local.id = 0; local.placed_uv = NULL;
        local.pose_theta = local.pose_x = local.pose_y = 0; one.charts = &local;
        for (size_t i = 0; i < opened.n_rels; i++) {
            const AsmRelation *r = opened.rels+i;
            if ((size_t)r->a == chart && r->a == r->b && (r->continuity & ASM_CONT_INTRACHART_CUT)) one.n_rels++;
        }
        one.rels = ARENA_ALLOC(scratch,one.n_rels*sizeof(AsmRelation)); one.corr = opened.corr; one.n_corr = opened.n_corr;
        size_t at = 0;
        for (size_t i = 0; i < opened.n_rels; i++) {
            const AsmRelation *r = opened.rels+i;
            if ((size_t)r->a == chart && r->a == r->b && (r->continuity & ASM_CONT_INTRACHART_CUT)) {
                one.rels[at] = *r; one.rels[at].a = one.rels[at].b = 0; at++;
            }
        }
        AsmField *f = AsmField_new_local(scratch,&one,NULL,NULL);
        if (!f) { stop = "intrinsic_cut_field_failed"; goto release; }
        memset(&result,0,sizeof result); result.before = amr_metric(f,0,f->uv);
        int refined = 0;
        if (budget.admission_solver_deadline > 0) {
#pragma omp critical(asm_region_deformation)
            { refined = amb_refine_region_budget(f,0,NULL,NULL,budget.deadline,precondition,&result,&budget,1); }
        } else refined = amb_refine_region_budget(f,0,NULL,NULL,budget.deadline,precondition,&result,&budget,1);
        result.after = amr_metric(f,0,f->uv);
        if (!refined) { stop = result.reason; goto release; }
        float *uv = ARENA_ALLOC(trial->arena,2*local.nv*sizeof(float)); double sign = local.flags & ASM_CHART_MIRROR ? -1 : 1;
        for (size_t v = 0; v < local.nv; v++) {
            uv[2*v] = (float)(sign*f->uv[2*v]); uv[2*v+1] = (float)f->uv[2*v+1];
            f->uv[2*v] = sign*(double)uv[2*v]; f->uv[2*v+1] = uv[2*v+1];
        }
        result.after = amr_metric(f,0,f->uv);
        if (!amr_qualified(&result.after)) { stop = "intrinsic_cut_float_metric_failed"; goto release; }
        if (!amb_cuts_pass(f)) { stop = "intrinsic_cut_gap_failed"; goto release; }
        if (apr_expired(&budget)) { stop = "intrinsic_cut_deadline"; goto release; }
        if (!f->contacts && AsmField_build_contacts(f)) { stop = "intrinsic_cut_contact_tree_failed"; goto release; }
        memset(&contacts,0,sizeof contacts);
        if (AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&contacts) || !contacts.complete || contacts.pairs) {
            stop = "intrinsic_cut_contact_failed"; goto release;
        }
        if (apr_expired(&budget) || !apr_workspace_check(&budget,0)) { stop = "intrinsic_cut_budget"; goto release; }
        opened.charts[chart].uv = uv; opened.charts[chart].placed_uv = NULL;
        opened.charts[chart].pose_theta = opened.charts[chart].pose_x = opened.charts[chart].pose_y = 0;
        *trial = opened; okay = 1; stop = tear ? "intrinsic_torn_cut_qualified" : "intrinsic_cut_qualified";
release:
        Arena_dispose(&scratch);
    }
done:
    *reason = stop;
    if (trace) fprintf(trace,"{\"phase\":\"intrinsic_cut_recovery\",\"chart\":%zu,\"reason\":\"%s\",\"qualified\":%s,"
        "\"cut_edges\":%zu,\"cut_length\":%.17g,\"before_bad\":%zu,\"after_bad\":%zu,\"steps\":%zu,"
        "\"cut_seconds\":%.17g,\"factor_seconds\":%.17g,\"seconds\":%.17g}\n",
        chart,stop,okay?"true":"false",cut.edges,cut.length,result.before.bad,result.after.bad,result.refine_steps,
        cut.seconds,result.factor_seconds,ves_clock_sec()-start);
    return okay;
}

static int api_cut_recover(AsmRun *trial, size_t chart, const AsmRepairOpts *opts, const char **reason, FILE *trace)
{ return api_cut_recover_from(trial,chart,opts,2,0,reason,trace,NULL); }

static void *api_copy(Arena_T arena, const void *source, size_t bytes)
{
    if (!source) return NULL; void *copy = ARENA_ALLOC(arena,bytes); memcpy(copy,source,bytes); return copy;
}

/* Prepare all ownership transfers before the overlay is published. Only
 * absent core charts may gain aliases; incumbent topology remains identical.
 * Existing measured observations keep every byte except their local alias
 * indices, which must still identify the same original source vertices. */
static int api_cut_stage(AsmRun *run, const AsmRun *trial, const uint8_t *core, AsmRun *staged)
{
    *staged = *run;
    if (trial->n_rels == run->n_rels && trial->n_corr == run->n_corr) {
        if (trial->n_charts != run->n_charts) return 0;
        for (size_t c = 0; c < run->n_charts; c++) {
            const AsmChart *a = run->charts+c, *b = trial->charts+c;
            if (a->nv != b->nv || a->nf != b->nf || a->xyz != b->xyz || a->faces != b->faces || a->vid != b->vid) return 0;
        }
        return 1;
    }
    if (trial->n_charts != run->n_charts || trial->n_rels < run->n_rels || trial->n_corr < run->n_corr ||
        trial->n_corr-run->n_corr != 2*(trial->n_rels-run->n_rels)) return 0;
    for (size_t c = 0; c < run->n_charts; c++) {
        const AsmChart *a = run->charts+c, *b = trial->charts+c;
        if (a->nv == b->nv && a->faces == b->faces) continue;
        if (!core[c] || AsmChart_registered(a) || b->nv < a->nv || a->id != b->id || a->comp != b->comp ||
            a->parent != b->parent || a->area3d != b->area3d || !acut_source_same(a,b) ||
            memcmp(a->xyz,b->xyz,3*a->nv*sizeof(float)) || memcmp(a->vid,b->vid,a->nv*sizeof(int32_t))) return 0;
    }
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels+i, *s = trial->rels+i;
        if (r->a != s->a || r->b != s->b || r->corr_first != s->corr_first || r->corr_count != s->corr_count) return 0;
        for (int k = 0; k < r->corr_count; k++) {
            const AsmCorr *a = run->corr+(size_t)r->corr_first+k, *b = trial->corr+(size_t)r->corr_first+k;
            AsmCorr restored = *b; restored.va = a->va; restored.vb = a->vb;
            if (memcmp(a,&restored,sizeof *a)) return 0;
            for (int side = 0; side < 2; side++) {
                size_t c = (size_t)(side?r->b:r->a); int32_t va = side?a->vb:a->va, vb = side?b->vb:b->va;
                if (va == vb) continue;
                if (!core[c] || va < 0 || vb < 0 || (size_t)va >= run->charts[c].nv || (size_t)vb >= trial->charts[c].nv ||
                    run->charts[c].vid[va] != trial->charts[c].vid[vb] ||
                    memcmp(run->charts[c].xyz+3*(size_t)va,trial->charts[c].xyz+3*(size_t)vb,3*sizeof(float))) return 0;
            }
        }
    }
    for (size_t i = run->n_rels; i < trial->n_rels; i++) {
        const AsmRelation *r = trial->rels+i;
        if (r->a < 0 || r->a != r->b || (size_t)r->a >= run->n_charts || !core[r->a] ||
            !(r->continuity & ASM_CONT_INTRACHART_CUT) || r->corr_first != (int32_t)(run->n_corr+2*(i-run->n_rels)) || r->corr_count != 2) return 0;
    }
    staged->charts = api_copy(run->arena,run->charts,run->n_charts*sizeof(AsmChart));
    for (size_t c = 0; c < run->n_charts; c++) if (trial->charts[c].nv != run->charts[c].nv) {
        const AsmChart *b = trial->charts+c; AsmChart *a = staged->charts+c; *a = *b;
        a->xyz = api_copy(run->arena,b->xyz,3*b->nv*sizeof(float)); a->nrm = api_copy(run->arena,b->nrm,3*b->nv*sizeof(float));
        a->vid = api_copy(run->arena,b->vid,b->nv*sizeof(int32_t)); a->faces = api_copy(run->arena,b->faces,3*b->nf*sizeof(int32_t));
        a->uv = api_copy(run->arena,b->uv,2*b->nv*sizeof(float)); a->placed_uv = NULL;
        a->boundary = api_copy(run->arena,b->boundary,b->nv); a->repaired = api_copy(run->arena,b->repaired,b->nv);
    }
    staged->n_rels = staged->cap_rels = trial->n_rels; staged->n_corr = staged->cap_corr = trial->n_corr;
    staged->rels = ARENA_ALLOC(run->arena,trial->n_rels*sizeof(AsmRelation));
    memcpy(staged->rels,run->rels,run->n_rels*sizeof(AsmRelation));
    memcpy(staged->rels+run->n_rels,trial->rels+run->n_rels,(trial->n_rels-run->n_rels)*sizeof(AsmRelation));
    staged->corr = api_copy(run->arena,trial->corr,trial->n_corr*sizeof(AsmCorr));
    return 1;
}
