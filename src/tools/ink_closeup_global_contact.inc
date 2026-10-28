/* Component feasibility against the complete current assembled material.
 * A refusal preserves the incumbent; this does not solve coupled contacts.
 * The query owns a separate arena: its lazy caches outlive CA scratch marks. */
typedef struct IcAlignGlobalContact {
    const MeshBinData *mesh;
    AsmContacts *query;
    double *incumbent, *trial;
    uint8_t *moving;
    const int32_t *vertices;
    size_t local_vertices, calls, refused, unknown;
    void *extra_context;
    int (*extra_guard)(void *, size_t, const double *);
} IcAlignGlobalContact;

static IcAlignGlobalContact *ic_align_global_new(Arena_T arena, const MeshBinData *mesh,
                                                double *incumbent, const int32_t *aliases)
{
    IcAlignGlobalContact *g = ARENA_CALLOC(arena,1,sizeof *g);
    AsmContactStats seed = {0};
    int threads = 1;
    g->mesh = mesh; g->incumbent = incumbent;
    g->trial = ARENA_ALLOC(arena,2*mesh->nv*sizeof *g->trial);
    memcpy(g->trial,incumbent,2*mesh->nv*sizeof *g->trial);
    g->moving = ARENA_CALLOC(arena,mesh->nf,sizeof *g->moving);
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    g->query = AsmContacts_new_source(arena,mesh->nv,mesh->faces,mesh->nf,g->trial,
                                     AsmContacts_material_tolerance(),threads,aliases);
    if (!g->query || AsmContacts_measure(g->query,g->trial,NULL,NULL,NULL,&seed) ||
        !ic_material_qualified(&seed)) {
        fprintf(stderr,"[global component guard] initial whole surface not certified: complete=%d pairs=%zu unknown=%zu\n",
                seed.complete,seed.pairs,seed.unmeasured);
        return NULL;
    }
    return g;
}

static int ic_align_global_candidate(void *context, size_t nv, const double *uv)
{
    IcAlignGlobalContact *g = context;
    if (!g || nv != g->local_vertices || !g->vertices || !uv) return 0;
    g->calls++;
    for (size_t k = 0; k < nv; k++) {
        size_t v = (size_t)g->vertices[k];
        g->trial[2*v] = uv[2*k]; g->trial[2*v+1] = uv[2*k+1];
    }
    AsmContactStats measured = {0};
    int rc = AsmContacts_measure(g->query,g->trial,NULL,NULL,NULL,&measured);
    int known = !rc && measured.complete && !measured.unmeasured &&
                isfinite(measured.area) && measured.area >= 0;
    int accept = known && !measured.pairs;
    /* The next call refits this same moving component. Before the mask changes,
     * commit() must also restore the cached boxes to the retained incumbent. */
    for (size_t k = 0; k < nv; k++) {
        size_t v = (size_t)g->vertices[k];
        g->trial[2*v] = g->incumbent[2*v]; g->trial[2*v+1] = g->incumbent[2*v+1];
    }
    if (accept && g->extra_guard) accept = g->extra_guard(g->extra_context,nv,uv) > 0;
    g->refused += !accept; g->unknown += !known;
    return accept;
}

/* The winner can differ from the last tested proposal, including when that
 * proposal was refused. Refit while the old component is still marked moving. */
static int ic_align_global_commit(IcAlignGlobalContact *g)
{
    if (!g || !g->query || !g->vertices) return -1;
    for (size_t k = 0; k < g->local_vertices; k++) {
        size_t v = (size_t)g->vertices[k];
        g->trial[2*v] = g->incumbent[2*v]; g->trial[2*v+1] = g->incumbent[2*v+1];
    }
    return AsmContacts_refit_changed(g->query,g->trial,g->moving);
}
