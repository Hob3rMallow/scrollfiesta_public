/* Strict placement certificates over original faces. A chart's triangle
 * cannot intersect the candidate outside their vertex bounding boxes, so
 * distant charts can be omitted without sampling their geometry. The index
 * uses authoritative global coordinates and is updated after every commit. */
typedef struct ApsBox { double lo[2], hi[2]; int have; } ApsBox;
struct AplStrict {
    ApsBox *boxes;
    size_t charts, candidate_cap, field_cap, checked, accepted;
    size_t walk_budget; /* zero retains the contact walk's normal finite cap */
    double deadline, seconds;
    FILE *log;
    int valid, io;
};
typedef struct ApsResult {
    size_t candidates, candidate_faces, field_charts, field_faces, bad_faces, bad_charts, failed_seams;
    AsmContactStats contacts;
    double seconds;
    const char *reason;
} ApsResult;

static int aps_box(const AsmChart *ch, ApsBox *box)
{
    memset(box, 0, sizeof *box);
    if (!AsmChart_registered(ch)) return 1;
    box->lo[0] = box->lo[1] = DBL_MAX; box->hi[0] = box->hi[1] = -DBL_MAX;
    for (size_t v = 0; v < ch->nv; v++) {
        double p[2]; AsmChart_point(ch, v, p);
        for (int d = 0; d < 2; d++) {
            if (!isfinite(p[d])) return 0;
            if (p[d] < box->lo[d]) box->lo[d] = p[d];
            if (p[d] > box->hi[d]) box->hi[d] = p[d];
        }
    }
    box->have = ch->nv != 0;
    return box->have;
}

static int aps_intersects(const ApsBox *a, const ApsBox *b)
{
    if (!a->have || !b->have) return 0;
    for (int d = 0; d < 2; d++) if (a->lo[d] > b->hi[d] + 1e-7 || b->lo[d] > a->hi[d] + 1e-7) return 0;
    return 1;
}

static AplStrict *aps_new(Arena_T arena, const AsmRun *run, const char *dir)
{
    AplStrict *s = ARENA_CALLOC(arena, 1, sizeof *s);
    s->charts = run->n_charts; s->candidate_cap = 300000; s->field_cap = 1000000;
    s->valid = s->io = 1; s->boxes = ARENA_CALLOC(arena, s->charts ? s->charts : 1, sizeof *s->boxes);
    for (size_t c = 0; c < s->charts; c++) if (!aps_box(run->charts + c, s->boxes + c)) s->valid = 0;
    if (dir) {
        char path[2048]; snprintf(path, sizeof path, "%s/stage5_strict.csv", dir); s->log = fopen(path, "wb");
        s->io = s->log && fprintf(s->log, "attempt,component,cluster,accepted,reason,candidate_charts,candidate_faces,field_charts,field_faces,bad_faces,bad_charts,failed_required_seams,overlap_pairs,contact_complete,seconds\n") >= 0;
    }
    return s;
}

static int aps_source_required(const AsmRelation *r, int selected)
{
    return !(r->flags & ASM_REL_CONTACT) && (r->continuity & (ASM_CONT_SOURCE | ASM_CONT_CUT)) &&
        (selected || AsmRel_is_join(r) || (r->continuity & (ASM_CONT_REQUIRED | ASM_CONT_CUT)));
}

static void aps_require_selected(AsmRun *run, const uint8_t *selected)
{
    for (size_t r = 0; r < run->n_rels; r++) if (selected[r] && aps_source_required(run->rels + r, 1))
        run->rels[r].continuity |= ASM_CONT_REQUIRED;
}

/* full is a selftest oracle: enumerate every represented face instead of
 * pruning chart boxes, keeping exactly the same seam and metric policy. */
static int aps_check(const AplCtx *x, int32_t component, int full, ApsResult *out)
{
    double started = ves_clock_sec(); Arena_T arena = Arena_new(); int okay = 0;
    const AplStrict *index = x->strict; size_t nc = x->run->n_charts;
    memset(out, 0, sizeof *out); out->reason = "strict_unavailable";
    if (!index || !index->valid || !index->io || index->charts != nc || !x->nproposal) goto done;
    if (index->deadline > 0 && ves_clock_sec() >= index->deadline) { out->reason = "strict_deadline"; goto done; }
    AsmRun trial = *x->run; trial.arena = arena;
    trial.charts = ARENA_ALLOC(arena, nc * sizeof(AsmChart)); memcpy(trial.charts, x->run->charts, nc * sizeof(AsmChart));
    uint8_t *core = ARENA_CALLOC(arena, nc, 1), *member = ARENA_CALLOC(arena, nc, 1);
    ApsBox bounds = {{DBL_MAX,DBL_MAX},{-DBL_MAX,-DBL_MAX},0};
    for (size_t k = 0; k < x->nproposal; k++) {
        const AsmChart *ch = x->proposal + k; int32_t c = ch->id;
        out->reason = "strict_invalid_proposal";
        if (c < 0 || (size_t)c >= nc || core[c] || x->orig_comp[c] != component || !AsmChart_registered(ch) || AsmChart_registered(x->run->charts + c)) goto done;
        trial.charts[c] = *ch; core[c] = member[c] = 1; out->candidates++; out->candidate_faces += ch->nf;
        ApsBox b; if (!aps_box(ch, &b)) goto done;
        for (int d = 0; d < 2; d++) {
            if (b.lo[d] < bounds.lo[d]) bounds.lo[d] = b.lo[d];
            if (b.hi[d] > bounds.hi[d]) bounds.hi[d] = b.hi[d];
        }
        bounds.have = 1;
    }
    for (size_t c = 0; c < nc; c++) if (x->orig_comp[c] == component && AsmChart_in_layout(trial.charts + c) && !core[c]) goto done;
    if (!full && out->candidate_faces > index->candidate_cap) { out->reason = "strict_candidate_cap"; goto done; }
    for (size_t c = 0; c < nc; c++) if (!core[c] && AsmChart_registered(trial.charts + c) &&
        (full || aps_intersects(&bounds, index->boxes + c))) member[c] = 1;
    for (size_t r = 0; r < trial.n_rels; r++) {
        const AsmRelation *rel = trial.rels + r;
        if (rel->a < 0 || rel->b < 0 || (size_t)rel->a >= nc || (size_t)rel->b >= nc) goto done;
        if ((!core[rel->a] && !core[rel->b]) || !aps_source_required(rel, x->seam_selected && x->seam_selected[r])) continue;
        if (!AsmChart_registered(trial.charts + rel->a) || !AsmChart_registered(trial.charts + rel->b)) continue;
        member[rel->a] = member[rel->b] = 1;
    }
    for (size_t c = 0; c < nc; c++) if (member[c]) { out->field_charts++; out->field_faces += trial.charts[c].nf; }
    if (!full && out->field_faces > index->field_cap) { out->reason = "strict_field_cap"; goto done; }
    AsmField *f = AsmField_new_local(arena, &trial, member, NULL);
    out->reason = "strict_field_failed"; if (!f) goto done;
    uint8_t *active = ARENA_CALLOC(arena, f->nf ? f->nf : 1, 1);
    for (size_t c = 0; c < nc; c++) if (core[c]) {
        double area = 0, good = 0; size_t bad = 0;
        if (f->face_offset[c] == SIZE_MAX) goto done;
        for (size_t k = 0; k < trial.charts[c].nf; k++) {
            size_t t = f->face_offset[c] + k; active[t] = 1; double j[4], lo, hi, det;
            AsmField_face(f, t, f->uv, j, &lo, &hi, &det);
            bad += !UvGuard_positive(f->uv, f->faces + 3*t) || !AsmMetric_within25(lo,hi);
            area += f->metric[t].area; if (AsmMetric_within10(lo,hi)) good += f->metric[t].area;
        }
        out->bad_faces += bad; out->bad_charts += bad || !(area > 0) || good < .95*area - 32*DBL_EPSILON*area;
    }
    if (out->bad_charts) { out->reason = "strict_original_metric"; goto done; }
    for (size_t i = 0; i < f->n_seams; i++) {
        const AsmFieldSeam *s = f->seams + i;
        if ((!core[s->a] && !core[s->b]) || !aps_source_required(trial.rels + s->relation, x->seam_selected && x->seam_selected[s->relation])) continue;
        if (!AsmChart_registered(trial.charts + s->a) || !AsmChart_registered(trial.charts + s->b)) continue;
        AsmSeamMeasure m; AsmField_seam(f, i, f->uv, &m); out->failed_seams += !m.pass;
    }
    if (out->failed_seams) { out->reason = "strict_required_seam"; goto done; }
    out->reason = "strict_contact_incomplete";
    if (AsmField_build_contacts(f)) goto done;
    AsmContacts_active_faces(f->contacts, active);
    if (index->walk_budget) AsmContacts_leaf_pair_budget(f->contacts, index->walk_budget);
    if (AsmContacts_measure(f->contacts, f->uv, NULL, NULL, NULL, &out->contacts)) goto done;
    if (out->contacts.pairs) { out->reason = "strict_overlap"; goto done; }
    if (index->deadline > 0 && ves_clock_sec() >= index->deadline) { out->reason = "strict_deadline"; goto done; }
    out->reason = "strict_qualified"; okay = 1;
done:
    out->seconds = ves_clock_sec() - started; Arena_dispose(&arena); return okay;
}

static int aps_candidate(AplCtx *x, int32_t component)
{
    ApsResult r; int okay = aps_check(x, component, 0, &r); AplStrict *s = x->strict;
    s->checked++; s->accepted += okay; s->seconds += r.seconds;
    if (s->log && (fprintf(s->log, "%d,%d,%d,%d,%s,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%d,%.9g\n",
        x->attempt, component, x->cluster, okay, r.reason, r.candidates, r.candidate_faces, r.field_charts, r.field_faces,
        r.bad_faces, r.bad_charts, r.failed_seams, r.contacts.pairs, r.contacts.complete, r.seconds) < 0 || fflush(s->log))) s->io = 0;
    if (!okay) x->last_reason = r.reason;
    return okay && s->io;
}
