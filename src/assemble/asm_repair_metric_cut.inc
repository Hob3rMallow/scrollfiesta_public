/* Open certified source obstructions along shortest original-edge paths
 * to an existing boundary. Face-corner components across UNCUT edges define
 * the duplicated vertices, preserving manifold fans and every original face.
 * Every opened edge receives a separately validated source-cut obligation. */
typedef struct AcutResult {
    size_t edges, added_vertices, before_obstructions, after_obstructions;
    double length, seconds, target_perimeter;
    size_t boundary_components;
    size_t orientation_before, orientation_after, workspace_estimate, sampled_bytes;
    int source_preserved;
    const char *reason;
} AcutResult;

static size_t acut_next(size_t corner)
{ return 3*(corner/3)+(corner+1)%3; }

static int32_t acut_slot(const int32_t *offset, const int32_t *target, int32_t a, int32_t b)
{
    int32_t lo = offset[a], hi = offset[(size_t)a+1];
    while (lo < hi) { int32_t mid = lo+(hi-lo)/2; if (target[mid] < b) lo = mid+1; else hi = mid; }
    return lo < offset[(size_t)a+1] && target[lo] == b ? lo : -1;
}

static double acut_length(const AsmChart *ch, size_t a, size_t b)
{
    double squared = 0;
    for (int d = 0; d < 3; d++) { double x = (double)ch->xyz[3*a+d]-ch->xyz[3*b+d]; squared += x*x; }
    return sqrt(squared);
}

static void *acut_alias(Arena_T arena, const void *source, size_t width, size_t count, const size_t *original)
{
    if (!source) return NULL;
    uint8_t *copy = ARENA_ALLOC(arena,width*count); const uint8_t *from = source;
    for (size_t v = 0; v < count; v++) memcpy(copy+width*v,from+width*original[v],width);
    return copy;
}

static int acut_source_same(const AsmChart *before, const AsmChart *after)
{
    if (before->nf != after->nf || before->cube != after->cube || !before->vid || !after->vid) return 0;
    for (size_t k = 0; k < 3*before->nf; k++) {
        int32_t a = before->faces[k], b = after->faces[k];
        if (a < 0 || b < 0 || (size_t)a >= before->nv || (size_t)b >= after->nv ||
            before->vid[a] != after->vid[b] || memcmp(before->xyz+3*(size_t)a,after->xyz+3*(size_t)b,3*sizeof(float))) return 0;
    }
    return 1;
}

static int acut_remap_point(const AsmChart *old, const AsmChart *next, const AsmTrim *trim, int32_t *vertex)
{
    if (!trim->valid) return 1;
    if (trim->face < 0 || (size_t)trim->face >= old->nf) return 0;
    for (int k = 0; k < 3; k++) if (old->faces[3*(size_t)trim->face+k] == *vertex) {
        *vertex = next->faces[3*(size_t)trim->face+k]; return 1;
    }
    return 0;
}

/* Cutting a boundary can change an inferred seam's arclength weights. Keep
 * the pre-cut obligations as the reference: a changed or incomplete old
 * witness is a refusal, not a newly weakened baseline for metric repair. */
static int acut_old_seams_budget(const AsmRun *run, const AsmRun *trial, size_t chart,
                                  size_t new_relations, const AsmRepairOpts *opts, size_t cut_bytes)
{
    Arena_T arena = Arena_new(); int okay = 0;
    if (run->n_charts != trial->n_charts || chart >= run->n_charts) goto done;
    uint8_t *member = ARENA_CALLOC(arena,run->n_charts,1); member[chart] = 1;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels+i;
        if (r->a < 0 || r->b < 0 || (size_t)r->a >= run->n_charts || (size_t)r->b >= run->n_charts) goto done;
        if ((size_t)r->a != chart && (size_t)r->b != chart) continue;
        member[r->a] = member[r->b] = 1;
    }
    size_t estimate = apr_member_workspace(trial,member);
    apr_bytes_add(&estimate,cut_bytes,1);
    if (apr_expired(opts) || !apr_workspace_check(opts,estimate)) goto done;
    AsmField *before = AsmField_new_local(arena,run,member,NULL), *after = AsmField_new_local(arena,trial,member,NULL);
    if (!before || !after || after->n_seams != before->n_seams+(AsmChart_registered(run->charts+chart)?new_relations:0) ||
        apr_expired(opts) || !apr_workspace_check(opts,estimate)) goto done;
    for (size_t i = 0; i < before->n_seams; i++) {
        const AsmFieldSeam *a = before->seams+i, *b = after->seams+i;
        if (a->relation != b->relation || a->count != b->count || a->complete != b->complete ||
            a->required != b->required || a->source_run != b->source_run || a->parity != b->parity || a->original_cut != b->original_cut) goto done;
        for (size_t k = 0; k < a->count; k++) {
            const AsmWitness *x = before->witness+a->first+k, *y = after->witness+b->first+k;
            if (x->valid != y->valid || x->weight != y->weight) goto done;
            if (x->valid) for (int side = 0; side < 2; side++) {
                double p[2], q[2]; AsmField_residual(x,side,before->uv,p); AsmField_residual(y,side,after->uv,q);
                if (memcmp(p,q,sizeof p)) goto done;
            }
        }
    }
    okay = 1;
done:
    Arena_dispose(&arena); return okay;
}

static int acut_old_seams_same(const AsmRun *run, const AsmRun *trial, size_t chart, size_t new_relations)
{
    AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1;
    return acut_old_seams_budget(run,trial,chart,new_relations,&opts,0);
}

static int acut_budget(const AsmRepairOpts *opts, size_t estimate, AcutResult *out)
{
    if (!apr_workspace_check(opts,estimate)) { out->reason = "workspace_memory_budget"; return 0; }
    if (apr_expired(opts)) { out->reason = "cut_deadline"; return 0; }
    return 1;
}

/* SIZE_MAX selects every proved positive-orientation obstruction. A
 * multi-source shortest-path forest has one existing-boundary root per tree.
 * Opening shared paths once keeps the chart connected, with no face removal. */
static int acut_make_mode_budget(const AsmRun *run, size_t chart, size_t vertex, Arena_T arena,
                                  AsmRun *trial, AcutResult *out, int outer, const AsmRepairOpts *opts)
{
    double start = ves_clock_sec(); memset(out,0,sizeof *out); out->reason = "invalid_chart"; int okay = 0;
    int automatic = vertex == SIZE_MAX;
    if (chart >= run->n_charts) goto done;
    const AsmChart *ch = run->charts+chart;
    if (outer < 0 || outer > 2 || (automatic && outer) || (!AsmChart_registered(ch) && !ch->uv) ||
        (!automatic && vertex >= ch->nv) || !ch->nv || !ch->nf || !ch->vid || !ch->xyz || !ch->faces) goto done;
    out->reason = "cut_index_range";
    if (ch->nv > INT32_MAX/2 || ch->nf > INT32_MAX/6 || run->n_rels > INT32_MAX || run->n_corr > INT32_MAX) goto done;
    /* All construction arrays, worst-case aliases and appended witnesses.
     * Nearby fields are preflighted separately once the actual forest exists. */
    size_t estimate = 0;
    apr_bytes_add(&estimate,ch->nv,160+2*(sizeof(AsmRelation)+2*sizeof(AsmCorr)));
    apr_bytes_add(&estimate,ch->nf,3*96+sizeof(AsmMetricFace));
    apr_bytes_add(&estimate,run->n_charts,sizeof(AsmChart));
    apr_bytes_add(&estimate,run->n_rels,sizeof(AsmRelation));
    apr_bytes_add(&estimate,run->n_corr,sizeof(AsmCorr));
    if (!acut_budget(opts,estimate,out)) goto done;
    uint8_t *obstruction = ARENA_CALLOC(arena,ch->nv,1); AmfStats before = amf_chart_marks(ch,NULL,obstruction);
    out->before_obstructions = before.obstruction_stars; out->orientation_before = before.orientation_obstructions;
    if (!before.io || (automatic ? !before.orientation_obstructions : !obstruction[vertex])) {
        out->reason = "vertex_has_no_proved_obstruction"; goto done;
    }
    if (!acut_budget(opts,estimate,out)) goto done;
    CSR_T graph = CSR_from_faces(arena,ch->faces,ch->nf,ch->nv);
    const int32_t *offset = CSR_offset(graph), *target = CSR_target(graph); size_t nnz = (size_t)CSR_nnz(graph);
    size_t *corner = ARENA_ALLOC(arena,nnz*sizeof(size_t));
    for (size_t k = 0; k < nnz; k++) corner[k] = SIZE_MAX;
    out->reason = "nonmanifold_source";
    for (size_t k = 0; k < 3*ch->nf; k++) {
        if (!(k%4096) && !acut_budget(opts,estimate,out)) goto done;
        int32_t slot = acut_slot(offset,target,ch->faces[k],ch->faces[acut_next(k)]);
        if (slot < 0 || corner[slot] != SIZE_MAX) goto done; corner[slot] = k;
    }
    uint8_t *boundary = ARENA_CALLOC(arena,ch->nv,1), *cut = ARENA_CALLOC(arena,nnz,1);
    size_t *boundary_root = ARENA_ALLOC(arena,ch->nv*sizeof(size_t));
    for (size_t v = 0; v < ch->nv; v++) boundary_root[v] = v;
    for (size_t v = 0; v < ch->nv; v++) for (int32_t k = offset[v]; k < offset[v+1]; k++) {
        int32_t reverse = acut_slot(offset,target,target[k],(int32_t)v); if (reverse < 0) goto done;
        if (corner[k] == SIZE_MAX || corner[reverse] == SIZE_MAX) {
            boundary[v] = boundary[target[k]] = 1;
            boundary_root[apr_root(boundary_root,v)] = apr_root(boundary_root,(size_t)target[k]);
        }
    }
    double *perimeter = ARENA_CALLOC(arena,ch->nv,sizeof(double)); size_t largest = SIZE_MAX;
    for (size_t v = 0; v < ch->nv; v++) for (int32_t k = offset[v]; k < offset[v+1]; k++)
        if (corner[k] == SIZE_MAX) perimeter[apr_root(boundary_root,v)] += acut_length(ch,v,(size_t)target[k]);
    for (size_t v = 0; v < ch->nv; v++) if (perimeter[v] > 0) {
        out->boundary_components++;
        if (largest == SIZE_MAX || perimeter[v] > perimeter[largest]) largest = v;
    }
    double *distance = ARENA_ALLOC(arena,ch->nv*sizeof(double));
    size_t *previous = ARENA_ALLOC(arena,ch->nv*sizeof(size_t)), *heap = ARENA_ALLOC(arena,ch->nv*sizeof(size_t));
    int32_t *position = ARENA_ALLOC(arena,ch->nv*sizeof(int32_t)); uint8_t *settled = ARENA_CALLOC(arena,ch->nv,1);
    for (size_t v = 0; v < ch->nv; v++) { distance[v] = DBL_MAX; previous[v] = SIZE_MAX; position[v] = -1; }
    size_t count = 0, end = SIZE_MAX, other = SIZE_MAX, remaining = before.orientation_obstructions, visits = 0;
    if (automatic) {
        for (size_t v = 0; v < ch->nv; v++) if (boundary[v]) {
            distance[v] = 0; position[v] = (int32_t)count; heap[count++] = v;
        }
    } else { count = 1; heap[0] = vertex; position[vertex] = 0; distance[vertex] = 0; }
    while (count) {
        if (!(visits++%4096) && !acut_budget(opts,estimate,out)) goto done;
        size_t v = acl_heap_pop(heap,position,&count,distance); settled[v] = 1;
        if (automatic && (obstruction[v] & AMF_POSITIVE_ORIENTATION) && !--remaining) break;
        if (!automatic && boundary[v]) {
            if (!outer) { end = v; break; }
            if (apr_root(boundary_root,v) == largest && end == SIZE_MAX) end = v;
            if (apr_root(boundary_root,v) != largest && other == SIZE_MAX) other = v;
            if (end != SIZE_MAX && (outer == 1 || other != SIZE_MAX)) break;
            continue; /* A cut may end at a boundary, never pass through it. */
        }
        for (int32_t k = offset[v]; k < offset[v+1]; k++) {
            size_t w = (size_t)target[k]; if (settled[w]) continue;
            if (outer == 1 && boundary[w] && apr_root(boundary_root,w) != largest) continue;
            double next = distance[v]+acut_length(ch,v,w); if (!(next < distance[w])) continue;
            distance[w] = next; previous[w] = v;
            if (position[w] < 0) { position[w] = (int32_t)count; heap[count++] = w; }
            acl_heap_up(heap,position,(size_t)position[w],distance);
        }
    }
    out->reason = "no_boundary_path";
    if (automatic ? remaining != 0 : (end == SIZE_MAX || end == vertex || (outer == 2 && other == SIZE_MAX))) goto done;
    int32_t *slots = ARENA_ALLOC(arena,ch->nv*sizeof(int32_t)); size_t path = 0;
    uint8_t *reached = ARENA_CALLOC(arena,ch->nv,1);
    for (size_t branch = 0; branch < (automatic ? ch->nv : outer == 2 ? 2 : 1); branch++) {
        if (automatic && !(obstruction[branch] & AMF_POSITIVE_ORIENTATION)) continue;
        size_t v = automatic ? branch : branch ? other : end;
        if (!automatic) reached[apr_root(boundary_root,v)] = 1;
        for (; automatic ? !boundary[v] : v != vertex; v = previous[v]) {
            if (previous[v] == SIZE_MAX) goto done;
            int32_t a = acut_slot(offset,target,(int32_t)previous[v],(int32_t)v), b = acut_slot(offset,target,(int32_t)v,(int32_t)previous[v]);
            if (a < 0 || b < 0 || corner[a] == SIZE_MAX || corner[b] == SIZE_MAX) goto done;
            if (cut[a]) break; /* Shared shortest-path prefix is opened once. */
            if (path >= ch->nv) { out->reason = "cut_forest_cycle"; goto done; }
            cut[a] = cut[b] = 1; slots[path++] = a; out->length += acut_length(ch,v,previous[v]);
        }
        if (automatic && boundary[v]) reached[apr_root(boundary_root,v)] = 1;
    }
    out->edges = path;
    for (size_t v = 0; v < ch->nv; v++) if (reached[v]) out->target_perimeter += perimeter[v];
    if (!path || path > INT32_MAX-run->n_rels || path > (INT32_MAX-run->n_corr)/2) { out->reason = "cut_index_range"; goto done; }
    if (!acut_budget(opts,estimate,out)) goto done;
    size_t ncorners = 3*ch->nf; size_t *parent = ARENA_ALLOC(arena,ncorners*sizeof(size_t));
    for (size_t k = 0; k < ncorners; k++) parent[k] = k;
    for (size_t v = 0; v < ch->nv; v++) for (int32_t k = offset[v]; k < offset[v+1]; k++) {
        if ((size_t)target[k] < v || cut[k] || corner[k] == SIZE_MAX) continue;
        int32_t reverse = acut_slot(offset,target,target[k],(int32_t)v);
        if (reverse < 0 || corner[reverse] == SIZE_MAX) continue;
        size_t a = corner[k], b = corner[reverse];
        parent[apr_root(parent,a)] = apr_root(parent,acut_next(b));
        parent[apr_root(parent,acut_next(a))] = apr_root(parent,b);
    }
    size_t *root_index = ARENA_ALLOC(arena,ncorners*sizeof(size_t));
    size_t *original = ARENA_ALLOC(arena,(ch->nv+ncorners)*sizeof(size_t));
    uint8_t *seen = ARENA_CALLOC(arena,ch->nv,1); size_t nv = ch->nv;
    for (size_t v = 0; v < ch->nv; v++) original[v] = v;
    for (size_t k = 0; k < ncorners; k++) root_index[k] = SIZE_MAX;
    int32_t *faces = ARENA_ALLOC(arena,ncorners*sizeof(int32_t));
    for (size_t k = 0; k < ncorners; k++) {
        if (!(k%4096) && !acut_budget(opts,estimate,out)) goto done;
        size_t root = apr_root(parent,k), old = (size_t)ch->faces[k];
        if (root_index[root] == SIZE_MAX) {
            size_t id = old;
            if (seen[old]) {
                if (nv >= INT32_MAX) { out->reason = "cut_index_range"; goto done; }
                id = nv++; original[id] = old;
            } else seen[old] = 1;
            root_index[root] = id;
        }
        faces[k] = (int32_t)root_index[root];
    }
    /* Joining two distinct boundary components adds one more alias than a
     * single tip-to-boundary slit. Extra aliases would split unrelated fans. */
    if (nv != ch->nv+path+(outer == 2)) { out->reason = "unexpected_cut_vertex_count"; goto done; }
    *trial = *run; trial->arena = arena; trial->cap_charts = run->n_charts;
    trial->charts = ARENA_ALLOC(arena,run->n_charts*sizeof(AsmChart)); memcpy(trial->charts,run->charts,run->n_charts*sizeof(AsmChart));
    AsmChart *changed = trial->charts+chart; changed->nv = nv; changed->faces = faces;
    changed->vid = acut_alias(arena,ch->vid,sizeof(int32_t),nv,original);
    changed->xyz = acut_alias(arena,ch->xyz,3*sizeof(float),nv,original);
    changed->nrm = acut_alias(arena,ch->nrm,3*sizeof(float),nv,original);
    changed->uv = acut_alias(arena,ch->uv,2*sizeof(float),nv,original);
    changed->repaired = acut_alias(arena,ch->repaired,1,nv,original);
    changed->boundary = acut_alias(arena,boundary,1,nv,original);
    changed->placed_uv = ARENA_ALLOC(arena,2*nv*sizeof(double));
    for (size_t v = 0; v < nv; v++) AsmChart_point(ch,original[v],changed->placed_uv+2*v);
    trial->n_rels = trial->cap_rels = run->n_rels+path; trial->n_corr = trial->cap_corr = run->n_corr+2*path;
    trial->rels = ARENA_CALLOC(arena,trial->n_rels,sizeof(AsmRelation));
    trial->corr = ARENA_CALLOC(arena,trial->n_corr,sizeof(AsmCorr));
    if (run->n_rels) memcpy(trial->rels,run->rels,run->n_rels*sizeof(AsmRelation));
    if (run->n_corr) memcpy(trial->corr,run->corr,run->n_corr*sizeof(AsmCorr));
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels+i;
        if ((size_t)r->a != chart && (size_t)r->b != chart) continue;
        for (int k = 0; k < r->corr_count; k++) {
            if (r->corr_first < 0 || (size_t)r->corr_first+(size_t)k >= run->n_corr) { out->reason = "invalid_existing_witness_range"; goto done; }
            AsmCorr *c = trial->corr+(size_t)r->corr_first+(size_t)k;
            if (((size_t)r->a == chart && !acut_remap_point(ch,changed,&c->trim_a,&c->va)) ||
                ((size_t)r->b == chart && !acut_remap_point(ch,changed,&c->trim_b,&c->vb))) { out->reason = "existing_witness_remap_failed"; goto done; }
        }
    }
    for (size_t k = 0; k < path; k++) {
        size_t left = corner[slots[k]]; int32_t reverse = acut_slot(offset,target,ch->faces[acut_next(left)],ch->faces[left]);
        size_t right = corner[reverse]; AsmRelation *r = trial->rels+run->n_rels+k;
        r->a = r->b = (int32_t)chart; r->flags = ASM_REL_WEAK;
        r->continuity = ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_REQUIRED | ASM_CONT_INTRACHART_CUT;
        r->corr_first = (int32_t)(run->n_corr+2*k); r->corr_count = 2; r->n_corr = 2;
        r->seam_len = acut_length(ch,(size_t)ch->faces[left],(size_t)ch->faces[acut_next(left)]);
        for (int side = 0; side < 2; side++) {
            AsmCorr *c = trial->corr+(size_t)r->corr_first+side;
            c->va = faces[side?acut_next(left):left]; c->vb = faces[side?right:acut_next(right)]; c->valid = 3;
            c->trim_a.valid = c->trim_b.valid = 1; c->trim_a.face = (int32_t)(left/3); c->trim_b.face = (int32_t)(right/3);
            changed->boundary[c->va] = changed->boundary[c->vb] = 1;
        }
    }
    out->added_vertices = nv-ch->nv; out->source_preserved = acut_source_same(ch,changed);
    AmfStats after = amf_chart(changed,NULL); out->after_obstructions = after.obstruction_stars;
    out->orientation_after = after.orientation_obstructions;
    if (!out->source_preserved || !after.io || after.obstruction_stars >= before.obstruction_stars) { out->reason = "cut_certificate_failed"; goto done; }
    if (automatic && after.orientation_obstructions) { out->reason = "orientation_obstruction_remaining"; goto done; }
    if (!acut_budget(opts,estimate,out)) goto done;
    if (!acut_old_seams_budget(run,trial,chart,path,opts,estimate)) {
        out->reason = "existing_seam_certificate_changed"; acut_budget(opts,estimate,out); goto done;
    }
    out->reason = "prepared"; okay = 1;
done:
    if (opts->workspace) { out->workspace_estimate = opts->workspace->estimate; out->sampled_bytes = opts->workspace->peak; }
    out->seconds = ves_clock_sec()-start; return okay;
}

static int acut_make_mode(const AsmRun *run, size_t chart, size_t vertex, Arena_T arena,
                           AsmRun *trial, AcutResult *out, int outer)
{
    AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1; opts.deadline = ves_clock_sec()+opts.patch_seconds;
    AprWorkspace work = {0}; work.measured = apr_process_bytes(&work.entry); opts.workspace = &work;
    return acut_make_mode_budget(run,chart,vertex,arena,trial,out,outer,&opts);
}

static int acut_make(const AsmRun *run, size_t chart, size_t vertex, Arena_T arena, AsmRun *trial, AcutResult *out)
{ return acut_make_mode(run,chart,vertex,arena,trial,out,0); }
