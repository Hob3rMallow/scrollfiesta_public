/* A known original edge does not need the six noisy correspondence samples
 * used to infer an unknown cross-cube seam. Its ENTIRE linear cut is measured
 * by its endpoints. Validate the source identity and both actual boundaries
 * first; the bit alone is never evidence. */
static int aff_cut_weights(const AsmChart *ch, const AsmRelation *r, const AsmCorr *corr, double weights[2])
{
    unsigned needed = ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_REQUIRED | ASM_CONT_INTRACHART_CUT;
    if (r->a != r->b || (r->continuity & needed) != needed || r->corr_count != 2 ||
        !ch->vid || !ch->xyz || !ch->faces || (r->flags & (ASM_REL_CONTACT | ASM_REL_PARITY | ASM_REL_DROPPED | ASM_REL_CROSSWRAP))) return 0;
    int32_t ids[2][2] = {{corr[0].va,corr[1].va},{corr[0].vb,corr[1].vb}};
    for (int k = 0; k < 2; k++) {
        if (corr[k].valid != 3 || ids[0][k] < 0 || ids[1][k] < 0 ||
            (size_t)ids[0][k] >= ch->nv || (size_t)ids[1][k] >= ch->nv) return 0;
        if (ch->vid[ids[0][k]] < 0 || ch->vid[ids[0][k]] != ch->vid[ids[1][k]] ||
            memcmp(ch->xyz+3*(size_t)ids[0][k],ch->xyz+3*(size_t)ids[1][k],3*sizeof(float))) return 0;
        for (int d = 0; d < 2; d++) if (corr[k].gap_a[d] != 0 || corr[k].gap_b[d] != 0) return 0;
    }
    if (ch->vid[ids[0][0]] == ch->vid[ids[0][1]] ||
        (ids[0][0] == ids[1][0] && ids[0][1] == ids[1][1])) return 0;
    int direction[2] = {0}, count[2] = {0}, original_count = 0;
    size_t face[2] = {SIZE_MAX,SIZE_MAX};
    for (size_t f = 0; f < ch->nf; f++) for (int k = 0; k < 3; k++) {
        int32_t a = ch->faces[3*f+k], b = ch->faces[3*f+(k+1)%3];
        if (a < 0 || b < 0 || (size_t)a >= ch->nv || (size_t)b >= ch->nv) return 0;
        int32_t oa = ch->vid[a], ob = ch->vid[b], first = ch->vid[ids[0][0]], second = ch->vid[ids[0][1]];
        original_count += (oa == first && ob == second) || (oa == second && ob == first);
        for (int side = 0; side < 2; side++) {
            int sign = a == ids[side][0] && b == ids[side][1] ? 1 : a == ids[side][1] && b == ids[side][0] ? -1 : 0;
            if (sign) { count[side]++; direction[side] = sign; face[side] = f; }
        }
    }
    if (original_count != 2 || count[0] != 1 || count[1] != 1 || direction[0] != -direction[1] || face[0] == face[1]) return 0;
    for (int k = 0; k < 2; k++) if (!corr[k].trim_a.valid || !corr[k].trim_b.valid ||
        corr[k].trim_a.face < 0 || corr[k].trim_b.face < 0 ||
        (size_t)corr[k].trim_a.face != face[0] || (size_t)corr[k].trim_b.face != face[1]) return 0;
    double squared = 0;
    for (int d = 0; d < 3; d++) { double v = (double)ch->xyz[3*(size_t)ids[0][1]+d]-ch->xyz[3*(size_t)ids[0][0]+d]; squared += v*v; }
    double length = sqrt(squared); if (!(length > 0) || !isfinite(length)) return 0;
    weights[0] = weights[1] = .5*length; return 1;
}

static int32_t aff_material_root(int32_t *parent, int32_t v)
{
    int32_t root = v; while (parent[root] != root) root = parent[root];
    while (parent[v] != v) { int32_t next = parent[v]; parent[v] = root; v = next; }
    return root;
}

/* Boundary quadrature belongs to the original material graph. Reconstruct
 * its vertex incidence from validated cuts before inferring ordinary seams;
 * their UV trim branches still refer to the actual opened face corners. */
static const AsmChart *aff_source_charts(AsmField *f)
{
    const AsmRun *run = f->run; int32_t *parent = NULL;
    uint8_t *changed = NULL;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels+i;
        if (!(r->continuity & ASM_CONT_INTRACHART_CUT) || r->a < 0 || r->a != r->b ||
            (size_t)r->a >= run->n_charts || f->vertex_offset[r->a] == SIZE_MAX ||
            r->corr_count != 2 || r->corr_first < 0 || (size_t)r->corr_first > run->n_corr ||
            2 > run->n_corr-(size_t)r->corr_first) continue;
        const AsmCorr *corr = run->corr+r->corr_first; double weights[2];
        if (!aff_cut_weights(run->charts+r->a,r,corr,weights)) continue;
        if (!parent) {
            parent = ARENA_ALLOC(f->arena,f->nv*sizeof(int32_t)); changed = ARENA_CALLOC(f->arena,run->n_charts,1);
            for (size_t v = 0; v < f->nv; v++) parent[v] = (int32_t)v;
        }
        changed[r->a] = 1; size_t offset = f->vertex_offset[r->a];
        for (int k = 0; k < 2; k++) {
            int32_t a = aff_material_root(parent,(int32_t)(offset+(size_t)corr[k].va));
            int32_t b = aff_material_root(parent,(int32_t)(offset+(size_t)corr[k].vb));
            parent[a>b?a:b] = a<b?a:b;
        }
    }
    if (!parent) return run->charts;
    for (size_t v = 0; v < f->nv; v++) parent[v] = aff_material_root(parent,(int32_t)v);
    f->material_vertex = parent;
    AsmChart *charts = ARENA_ALLOC(f->arena,run->n_charts*sizeof(AsmChart));
    memcpy(charts,run->charts,run->n_charts*sizeof(AsmChart));
    for (size_t c = 0; c < run->n_charts; c++) if (changed[c]) {
        AsmChart *ch = charts+c; size_t offset = f->vertex_offset[c];
        const int32_t *faces = ch->faces; ch->faces = ARENA_ALLOC(f->arena,3*ch->nf*sizeof(int32_t));
        for (size_t k = 0; k < 3*ch->nf; k++) ch->faces[k] = parent[offset+(size_t)faces[k]]-(int32_t)offset;
    }
    return charts;
}

/* `loops` caches each source chart's boundary loops (in `cache`) across the
 * seams of one field build; the weights are identical to building them per seam. */
static int aff_boundary_weights(AsmField *f, const AsmChart *source, Arena_T cache, AsmBoundaryLoops **loops,
                                 const AsmRelation *r, const AsmCorr *corr, size_t count, double *weights)
{
    if (!loops[r->a]) loops[r->a] = AsmBoundary_loops(cache,source+r->a);
    if (!loops[r->b]) loops[r->b] = AsmBoundary_loops(cache,source+r->b);
    const AsmBoundaryLoops *la = loops[r->a], *lb = loops[r->b];
    if (!f->material_vertex) return AsmBoundary_weights_loops(f->arena,source+r->a,la,source+r->b,lb,corr,count,weights);
    Arena_Mark mark = Arena_save(f->arena);
    AsmCorr *mapped = ARENA_ALLOC(f->arena,count*sizeof(AsmCorr)); memcpy(mapped,corr,count*sizeof(AsmCorr));
    size_t a = f->vertex_offset[r->a], b = f->vertex_offset[r->b]; int okay = -1;
    for (size_t k = 0; k < count; k++) {
        if (corr[k].va < 0 || corr[k].vb < 0 || (size_t)corr[k].va >= source[r->a].nv || (size_t)corr[k].vb >= source[r->b].nv) goto done;
        mapped[k].va = f->material_vertex[a+(size_t)corr[k].va]-(int32_t)a;
        mapped[k].vb = f->material_vertex[b+(size_t)corr[k].vb]-(int32_t)b;
    }
    okay = AsmBoundary_weights_loops(f->arena,source+r->a,la,source+r->b,lb,mapped,count,weights);
done:
    Arena_restore(f->arena,mark); return okay;
}

/* Only independently validated cut witnesses may restore material adjacency.
 * Equal source IDs alone are insufficient: unrelated charts can reuse IDs.
 * The no-cut path allocates nothing and keeps the existing boundary index.
 * A torn cut (ASM_REL_TORN) is no seam obligation, but its material stays
 * continuous: aff_source_charts joins every validated intra-chart cut, and
 * that map is used even when no cut seam remains. A field holding only torn
 * cuts otherwise treated the tear as an original boundary and exempted a fold
 * across it, which the whole-run audit counted (PHerc0343_c chart 514: 0
 * pairs in its leaf field, 84 pairs and 46.9 vox^2 in the audit). */
static int aff_material_vertices(AsmField *f, int32_t **result)
{
    *result = NULL; size_t cuts = 0;
    for (size_t i = 0; i < f->n_seams; i++) {
        const AsmFieldSeam *s = f->seams+i;
        if (!s->original_cut || f->vertex_offset[s->a] == SIZE_MAX) continue;
        if (!s->complete || !s->parity || !s->source_run || !s->required || s->a != s->b || s->count != 2) return 0;
        for (size_t k = 0; k < s->count; k++) if (!f->witness[s->first+k].valid) return 0;
        cuts++;
    }
    if (f->material_vertex) { *result = f->material_vertex; return 1; }
    if (!cuts) return 1;
    int32_t *parent = ARENA_ALLOC(f->arena,f->nv*sizeof(int32_t));
    for (size_t v = 0; v < f->nv; v++) parent[v] = (int32_t)v;
    for (size_t i = 0; i < f->n_seams; i++) {
        const AsmFieldSeam *s = f->seams+i;
        if (!s->original_cut || f->vertex_offset[s->a] == SIZE_MAX) continue;
        const AsmRelation *r = f->run->rels+s->relation;
        for (size_t k = 0; k < 2; k++) {
            const AsmCorr *c = f->run->corr+(size_t)r->corr_first+k;
            int32_t a = (int32_t)(f->vertex_offset[s->a]+(size_t)c->va), b = (int32_t)(f->vertex_offset[s->b]+(size_t)c->vb);
            a = aff_material_root(parent,a); b = aff_material_root(parent,b);
            parent[a>b?a:b] = a<b?a:b;
        }
    }
    for (size_t v = 0; v < f->nv; v++) parent[v] = aff_material_root(parent,(int32_t)v);
    *result = parent; return 1;
}
