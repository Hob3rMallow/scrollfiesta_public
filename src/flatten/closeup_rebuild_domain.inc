/* Re-open the bounded continuation region erased by earlier failed inference.
 * This computes a proposal domain, not permission to join physical sheets.
 * Rebuild still checks actual source groups and DATA witnesses at every cell. */
int CloseupSurface_rebuild_regrowth_domain(const CloseupGrid *g, uint8_t *domain,
                                           double max_distance_px, CloseupRebuildDomainStats *st)
{
    uint8_t *seed = NULL; float *distance2 = NULL; int32_t *nearest = NULL;
    size_t nn = 0; int rc = -1;
    if (!g || !domain || !st || g->ni < 3 || g->nj < 3 || g->ni > INT32_MAX/g->nj ||
        !(g->h > 0) || !isfinite(g->h) || !(max_distance_px > 0) || !isfinite(max_distance_px)) return -1;
    memset(st,0,sizeof *st); nn = g->ni*g->nj;
    seed = (uint8_t *)calloc(nn,1); distance2 = (float *)malloc(nn*sizeof(float));
    nearest = (int32_t *)malloc(nn*sizeof(int32_t));
    if (!seed || !distance2 || !nearest) goto done;
    st->reach_px = max_distance_px;
    st->source_frame_bounds[0] = g->u0+g->h; st->source_frame_bounds[1] = g->v0+g->h;
    st->source_frame_bounds[2] = g->u0+g->h*(double)(g->ni-2);
    st->source_frame_bounds[3] = g->v0+g->h*(double)(g->nj-2);
    for (size_t k = 0; k < nn; k++) {
        size_t i = k%g->ni, j = k/g->ni;
        if (domain[k] && (!i || !j || i+1 == g->ni || j+1 == g->nj)) goto done;
        if (g->state[k] && !domain[k]) goto done;
        st->required_nodes += domain[k] != 0;
        st->required_enclosed_gap_nodes += domain[k] && !g->state[k];
        seed[k] = domain[k] && g->state[k] == CLOSEUP_NODE_DATA;
        st->data_donors += seed[k] != 0;
    }
    if (cs_edt2d(seed,g->ni,g->nj,distance2,nearest)) goto done;
    for (size_t k = 0; k < nn; k++) {
        size_t i = k%g->ni, j = k/g->ni;
        int reached = nearest[k] >= 0 && g->h*sqrt((double)distance2[k]) <= max_distance_px;
        if (domain[k]) st->required_beyond_reach += !reached;
        else if (i && j && i+1 < g->ni && j+1 < g->nj && reached) { domain[k] = 1; st->added_nodes++; }
        st->nodes += domain[k] != 0;
    }
    rc = 0;
done:
    free(seed); free(distance2); free(nearest); return rc;
}

static int cs_rebuild_regrowth_domain_selftest(void)
{
    CloseupGrid g = {0}; CloseupRebuildDomainStats st = {0}; uint8_t domain[49] = {0};
    int fails = 0;
    if (CloseupGrid_alloc(&g,0,12,0,12,2)) return 1;
    /* Active DATA at (2,2); a historical required island at (10,10) is kept
     * even beyond reach. A synthetic unused point is never imported as DATA. */
    domain[8] = domain[40] = 1;
    g.state[8] = CLOSEUP_NODE_DATA; g.state[40] = CLOSEUP_NODE_FILL;
    if (CloseupSurface_rebuild_regrowth_domain(&g,domain,4,&st)) fails++;
    if (!domain[10] || domain[11] || !domain[40] || st.required_nodes != 2 || st.required_beyond_reach != 1 ||
        st.data_donors != 1 || st.nodes != st.required_nodes+st.added_nodes) fails++;
    for (size_t k = 0; k < 49; k++) if ((k%g.ni == 0 || k%g.ni == g.ni-1 || k/g.ni == 0 || k/g.ni == g.nj-1) && domain[k]) fails++;
    if (g.state[8] != CLOSEUP_NODE_DATA || g.state[40] != CLOSEUP_NODE_FILL) fails++;
    memset(domain,0,sizeof domain); memset(g.state,0,49); domain[40] = 1; g.state[40] = CLOSEUP_NODE_FILL;
    if (CloseupSurface_rebuild_regrowth_domain(&g,domain,4,&st) || st.data_donors || st.added_nodes || st.nodes != 1) fails++;
    fprintf(stderr,"  DATA rebuild bounded regrowth domain / exact required union / DATA-only donors: %s\n",fails ? "FAIL" : "ok");
    CloseupGrid_free(&g); return fails != 0;
}
