/* Repair a completed regular-grid mesh without resampling retained material.
 * The first SLIM pass supplies the actual metric to Cartesian hole insertion;
 * the second flattens the enlarged surface. Original vertices and faces are
 * immutable prefixes. Cuts may duplicate them, with composed provenance. */
typedef struct IcFixedCompletion {
    MeshBinData mesh;
    float *nrm;
    uint8_t *attr;
    int32_t *chart, *source, *face_patch, *material;
    size_t npatches;
    double *frame;
} IcFixedCompletion;

static void ic_fixed_free(IcFixedCompletion *r)
{
    MeshBin_dispose(&r->mesh);
    free(r->nrm); free(r->attr); free(r->chart); free(r->source); free(r->frame); free(r->face_patch);
    free(r->material);
    memset(r,0,sizeof *r);
}

static int32_t ic_fixed_root(int32_t *parent, int32_t v)
{
    while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; }
    return v;
}

/* Triangle types in a cell: abc, abd, adc, bdc. The two diagonal pairs
 * are {abc,bdc} and {abd,adc}; a retained three-corner face owns its diagonal. */
static int ic_fixed_type(unsigned corners)
{
    return corners == 7 ? 1 : corners == 11 ? 2 : corners == 13 ? 4 : corners == 14 ? 8 : 0;
}

static size_t ic_fixed_faces(const CloseupGrid *g, const double *uv, const int32_t *vertex,
                             const uint8_t *old_cells, const int8_t *sign,
                             size_t original_nv, int32_t *faces)
{
    static const int triangle[4][3] = {{0,1,2},{0,1,3},{0,3,2},{1,3,2}};
    size_t count = 0;
    for (size_t j = 0; j+1 < g->nj; j++) for (size_t i = 0; i+1 < g->ni; i++) {
        size_t a = j*g->ni+i, nodes[4] = {a,a+1,a+g->ni,a+g->ni+1};
        unsigned old = old_cells[j*(g->ni-1)+i], valid = 0, fresh = 0, use = 0;
        if (old == 9 || old == 6) continue;
        for (int k = 0; k < 4; k++) if (g->state[nodes[k]]) {
            valid |= 1u<<k;
            if (vertex[nodes[k]] >= 0 && (size_t)vertex[nodes[k]] >= original_nv) fresh |= 1u<<k;
        }
        if (!fresh) continue;
        if (old) use = old == 1 ? 8 : old == 8 ? 1 : old == 2 ? 4 : old == 4 ? 2 : 0;
        else if (valid == 15) use = fresh == 1 || fresh == 8 ? 9 : 6;
        else use = (unsigned)ic_fixed_type(valid);
        for (int t = 0; t < 4; t++) if (use & (1u<<t)) {
            size_t v[3]; int touches = 0, ready = 1;
            for (int k = 0; k < 3; k++) {
                v[k] = nodes[triangle[t][k]];
                ready &= g->state[v[k]] && vertex[v[k]] >= 0;
                touches |= vertex[v[k]] >= 0 && (size_t)vertex[v[k]] >= original_nv;
            }
            if (!ready || !touches) continue;
            int32_t component = g->component[v[0]];
            if (component < 0 || component != g->component[v[1]] || component != g->component[v[2]]) continue;
            double det = (uv[2*v[1]]-uv[2*v[0]])*(uv[2*v[2]+1]-uv[2*v[0]+1])-
                         (uv[2*v[1]+1]-uv[2*v[0]+1])*(uv[2*v[2]]-uv[2*v[0]]);
            if (!(det > 0) || !isfinite(det) || !sign[component]) return SIZE_MAX;
            if (faces) {
                faces[3*count] = vertex[v[0]];
                faces[3*count+1] = vertex[v[sign[component] > 0 ? 1 : 2]];
                faces[3*count+2] = vertex[v[sign[component] > 0 ? 2 : 1]];
            }
            count++;
        }
    }
    return count;
}

static int ic_fixed_append(const MeshBinData *m, const double *pre_uv, const float *nrm,
                            const uint8_t *attr, const int32_t *chart,
                            IcFixedCompletion *r, CloseupFitStats *stats, CloseupProposalStats *proposal)
{
    CloseupGrid g = {0};
    int32_t *parent = NULL, *node = NULL, *vertex = NULL, *patch_nodes = NULL, *patch_vertices = NULL;
    int8_t *sign = NULL;
    uint8_t *active = NULL, *cells = NULL;
    double *grid_uv = NULL;
    double lo[2] = {INFINITY,INFINITY}, hi[2] = {-INFINITY,-INFINITY};
    size_t nn = 0, added = 0, added_faces = 0;
    const char *phase = "source topology";
    int rc = -1;
    if (!m->nv || !m->nf || m->nv > INT32_MAX || (!m->uv && !m->uv64)) goto done;
    parent = malloc(m->nv*sizeof *parent); node = malloc(m->nv*sizeof *node);
    sign = calloc(m->nv,sizeof *sign); active = calloc(m->nv,1);
    if (!parent || !node || !sign || !active) goto done;
    for (size_t v = 0; v < m->nv; v++) { parent[v] = (int32_t)v; node[v] = -1; }
    for (size_t f = 0; f < m->nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t v = m->faces[3*f+(size_t)k];
            if (v < 0 || (size_t)v >= m->nv) goto done;
            active[v] = 1;
        }
        for (int k = 1; k < 3; k++) {
            int32_t a = ic_fixed_root(parent,m->faces[3*f]), b = ic_fixed_root(parent,m->faces[3*f+(size_t)k]);
            if (a != b) parent[b] = a;
        }
    }
    for (size_t v = 0; v < m->nv; v++) if (active[v]) {
        parent[v] = ic_fixed_root(parent,(int32_t)v);
        for (size_t k = 0; k < 2; k++) {
            double p = MeshBin_uv(m,2*v+k);
            if (!isfinite(p) || !isfinite(pre_uv[2*v+k])) goto done;
            lo[k] = fmin(lo[k],p); hi[k] = fmax(hi[k],p);
        }
    }
    phase = "regular-grid import";
    if (CloseupGrid_alloc(&g,0,hi[0]-lo[0]+2*CLOSEUP_GRID_PX,
                             0,hi[1]-lo[1]+2*CLOSEUP_GRID_PX,CLOSEUP_GRID_PX)) goto done;
    g.u0 = lo[0]-g.h; g.v0 = lo[1]-g.h;
    if (g.ni > INT32_MAX/g.nj) goto done;
    nn = g.ni*g.nj;
    vertex = malloc(nn*sizeof *vertex); grid_uv = malloc(2*nn*sizeof *grid_uv);
    cells = calloc((g.ni-1)*(g.nj-1),1);
    if (!vertex || !grid_uv || !cells) goto done;
    for (size_t k = 0; k < nn; k++) { vertex[k] = -1; grid_uv[2*k] = grid_uv[2*k+1] = NAN; }
    for (size_t v = 0; v < m->nv; v++) if (active[v]) {
        double x = (MeshBin_uv(m,2*v)-g.u0)/g.h, y = (MeshBin_uv(m,2*v+1)-g.v0)/g.h;
        double ix = floor(x+.5), iy = floor(y+.5);
        if (fabs(x-ix) > 1e-5 || fabs(y-iy) > 1e-5 || ix < 0 || iy < 0 ||
            ix >= (double)g.ni || iy >= (double)g.nj) goto done;
        size_t at = (size_t)iy*g.ni+(size_t)ix;
        if (vertex[at] >= 0 || (attr[v] != 1 && attr[v] != 2 && attr[v] != 3)) goto done;
        node[v] = (int32_t)at; vertex[at] = (int32_t)v;
        memcpy(g.pos+3*at,m->verts+3*v,3*sizeof(float)); memcpy(g.nrm+3*at,nrm+3*v,3*sizeof(float));
        g.state[at] = attr[v] == 2 ? CLOSEUP_NODE_FILL : CLOSEUP_NODE_DATA;
        g.chart[at] = chart[v]; g.component[at] = parent[v];
        grid_uv[2*at] = pre_uv[2*v]; grid_uv[2*at+1] = pre_uv[2*v+1];
    }
    phase = "retained cell diagonals";
    for (size_t f = 0; f < m->nf; f++) {
        size_t x[3], y[3], x0 = g.ni, y0 = g.nj;
        for (int k = 0; k < 3; k++) {
            size_t at = (size_t)node[m->faces[3*f+(size_t)k]];
            x[k] = at%g.ni; y[k] = at/g.ni;
            if (x[k] < x0) x0 = x[k]; if (y[k] < y0) y0 = y[k];
        }
        unsigned corners = 0;
        for (int k = 0; k < 3; k++) {
            if (x[k]-x0 > 1 || y[k]-y0 > 1) goto done;
            corners |= 1u << (x[k]-x0+2*(y[k]-y0));
        }
        unsigned type = (unsigned)ic_fixed_type(corners);
        if (!type || x0+1 >= g.ni || y0+1 >= g.nj) goto done;
        uint8_t *cell = cells+y0*(g.ni-1)+x0;
        if (*cell & type) goto done;
        *cell |= (uint8_t)type;
        if (*cell != 1 && *cell != 2 && *cell != 4 && *cell != 8 && *cell != 6 && *cell != 9) goto done;
        double det = ((double)x[1]-(double)x[0])*((double)y[2]-(double)y[0])-
                     ((double)y[1]-(double)y[0])*((double)x[2]-(double)x[0]);
        int8_t s = det > 0 ? 1 : -1;
        int32_t component = parent[m->faces[3*f]];
        if (sign[component] && sign[component] != s) goto done;
        sign[component] = s;
    }
    phase = "fixed Cartesian completion";
    if (proposal) {
        if (CloseupSurface_propose_fixed(&g,CLOSEUP_COMPLETE_FILL_PX,&patch_nodes,proposal)) goto done;
        stats->variational_holes = proposal->holes; stats->variational_filled = proposal->proposed;
        stats->variational_refused = proposal->refused; stats->variational_nodes = proposal->nodes;
        r->npatches = proposal->proposed;
        for (size_t k = 0; k < nn; k++) if (g.state[k] && vertex[k] < 0) {
            grid_uv[2*k] = g.u0+g.h*(double)(k%g.ni);
            grid_uv[2*k+1] = g.v0+g.h*(double)(k/g.ni);
        }
    } else if (CloseupSurface_complete_fixed_uv(&g,grid_uv,CLOSEUP_COMPLETE_FILL_PX,stats)) goto done;
    for (size_t k = 0; k < nn; k++) if (g.state[k] && vertex[k] < 0) {
        if (m->nv+added >= INT32_MAX) goto done;
        vertex[k] = (int32_t)(m->nv+added++);
    }
    added_faces = ic_fixed_faces(&g,grid_uv,vertex,cells,sign,m->nv,NULL);
    if (added_faces == SIZE_MAX || m->nf+added_faces > SIZE_MAX/(3*sizeof(int32_t))) goto done;
    r->mesh.nv = m->nv+added; r->mesh.nf = m->nf+added_faces;
    r->mesh.verts = malloc(3*r->mesh.nv*sizeof(float)); r->mesh.uv64 = malloc(2*r->mesh.nv*sizeof(double));
    r->mesh.faces = malloc(3*r->mesh.nf*sizeof(int32_t)); r->nrm = malloc(3*r->mesh.nv*sizeof(float));
    r->attr = malloc(r->mesh.nv); r->chart = malloc(r->mesh.nv*sizeof(int32_t));
    r->source = malloc(r->mesh.nv*sizeof(int32_t)); r->frame = malloc(2*r->mesh.nv*sizeof(double));
    if (!r->mesh.verts || !r->mesh.uv64 || !r->mesh.faces || !r->nrm || !r->attr || !r->chart || !r->source || !r->frame) goto done;
    if (proposal) {
        r->face_patch = malloc(r->mesh.nf*sizeof(int32_t)); patch_vertices = malloc(r->mesh.nv*sizeof(int32_t));
        if (!r->face_patch || !patch_vertices) goto done;
        for (size_t f = 0; f < r->mesh.nf; f++) r->face_patch[f] = -1;
        for (size_t v = 0; v < r->mesh.nv; v++) patch_vertices[v] = -1;
    }
    memcpy(r->mesh.verts,m->verts,3*m->nv*sizeof(float)); memcpy(r->mesh.faces,m->faces,3*m->nf*sizeof(int32_t));
    memcpy(r->mesh.uv64,pre_uv,2*m->nv*sizeof(double)); memcpy(r->nrm,nrm,3*m->nv*sizeof(float));
    memcpy(r->attr,attr,m->nv); memcpy(r->chart,chart,m->nv*sizeof(int32_t));
    for (size_t v = 0; v < m->nv; v++) {
        r->source[v] = (int32_t)v;
        r->frame[2*v] = MeshBin_uv(m,2*v); r->frame[2*v+1] = MeshBin_uv(m,2*v+1);
    }
    for (size_t k = 0; k < nn; k++) if (vertex[k] >= 0 && (size_t)vertex[k] >= m->nv) {
        size_t v = (size_t)vertex[k];
        memcpy(r->mesh.verts+3*v,g.pos+3*k,3*sizeof(float)); memcpy(r->nrm+3*v,g.nrm+3*k,3*sizeof(float));
        memcpy(r->mesh.uv64+2*v,grid_uv+2*k,2*sizeof(double));
        r->attr[v] = CLOSEUP_NODE_FILL; r->chart[v] = -1; r->source[v] = -1;
        if (proposal) patch_vertices[v] = patch_nodes[k];
        r->frame[2*v] = g.u0+g.h*(double)(k%g.ni); r->frame[2*v+1] = g.v0+g.h*(double)(k/g.ni);
    }
    if (ic_fixed_faces(&g,grid_uv,vertex,cells,sign,m->nv,r->mesh.faces+3*m->nf) != added_faces) goto done;
    if (proposal) for (size_t f = m->nf; f < r->mesh.nf; f++) {
        int32_t p = -1;
        for (int k = 0; k < 3; k++) {
            int32_t q = patch_vertices[r->mesh.faces[3*f+(size_t)k]];
            if (q < 0) continue;
            if (p >= 0 && p != q) goto done;
            p = q;
        }
        if (p < 0 || (size_t)p >= r->npatches) goto done;
        r->face_patch[f] = p;
    }
    fprintf(stderr,"[ink_closeup] fixed completion: %zu/%zu holes filled, %zu refused; appended %zu vertices / %zu faces\n",
            stats->variational_filled,stats->variational_holes,stats->variational_refused,added,added_faces);
    rc = 0;
done:
    if (rc) { fprintf(stderr,"ERROR: retained-mesh completion failed at %s\n",phase); ic_fixed_free(r); }
    free(parent); free(node); free(sign); free(active); free(vertex); free(cells); free(grid_uv);
    free(patch_nodes); free(patch_vertices); CloseupGrid_free(&g);
    return rc;
}

static int ic_fixed_cut(IcFixedCompletion *r, const SlimRefineDirection *directions,
                         size_t *cones, size_t *edges, size_t *fiber_cones)
{
    int32_t *map = NULL;
    if (ic_cone_cut_directed(&r->mesh,directions,&map,cones,edges,fiber_cones)) return -1;
    if (!map) return 0;
    int32_t *material = malloc(r->mesh.nv*sizeof *material);
    if (!material) { free(map); return -1; }
    for (size_t v = 0; v < r->mesh.nv; v++) material[v] = r->material ? r->material[map[v]] : map[v];
    free(r->material); r->material = material;
    size_t stride[5] = {3*sizeof(float),1,sizeof(int32_t),sizeof(int32_t),2*sizeof(double)};
    void *old[5] = {r->nrm,r->attr,r->chart,r->source,r->frame}, *next[5] = {0};
    for (int k = 0; k < 5; k++) {
        next[k] = malloc(r->mesh.nv*stride[k]);
        if (!next[k]) { for (int j = 0; j < 5; j++) free(next[j]); free(map); return -1; }
        for (size_t v = 0; v < r->mesh.nv; v++)
            memcpy((uint8_t *)next[k]+v*stride[k],(uint8_t *)old[k]+(size_t)map[v]*stride[k],stride[k]);
    }
    for (int k = 0; k < 5; k++) free(old[k]);
    r->nrm = next[0]; r->attr = next[1]; r->chart = next[2]; r->source = next[3]; r->frame = next[4];
    free(map); return 0;
}

static int ic_stage_repair_completed(const IcConfig *c, const char *prepared_uv, int provisional)
{
    MeshBinData source = {0}, prepared = {0}; IcFixedCompletion result = {0};
    SlimRefineStats pre = {0}, final = {0}; CloseupFitStats fill = {0};
    AsmContactStats contact = {0};
    CloseupFiberStats measured = {0}; SlimRefineDirection *directions = NULL;
    CloseupProposalStats proposal = {0}; CloseupProposalMetricStats metric = {0};
    uint8_t *patch_pass = NULL;
    char path[2600], target[2600], json_source[4800];
    const char *suffix[3] = {"_surface_nrm.f32","_surface_attr.u8","_surface_chart.i32"};
    void *side[3] = {0}; size_t stride[3] = {3*sizeof(float),1,sizeof(int32_t)}, got = 0;
    double *uv = NULL; size_t cones = 0, cut_edges = 0, fiber_cones = 0, appended_nv = 0, appended_nf = 0;
    int rc = 1, mode = ic_reparameterize_mode; FILE *report = NULL;
    snprintf(path,sizeof path,"%s_surface.vmesh",c->surface_source);
    if (MeshBin_read_precise_malloc(path,&source)) goto done;
    for (int k = 0; k < 3; k++) {
        snprintf(path,sizeof path,"%s%s",c->surface_source,suffix[k]);
        side[k] = ic_read_file(path,source.nv*stride[k],&got);
        if (!side[k] || got != source.nv*stride[k]) goto done;
    }
    uv = malloc(2*source.nv*sizeof(double)); if (!uv) goto done;
    if (provisional) {
        if (prepared_uv) goto done;
        for (size_t k = 0; k < 2*source.nv; k++) uv[k] = MeshBin_uv(&source,k);
    } else if (prepared_uv) {
        if (MeshBin_read_precise_malloc(prepared_uv,&prepared) || prepared.nv != source.nv ||
            prepared.nf != source.nf || (!prepared.uv && !prepared.uv64) ||
            memcmp(prepared.verts,source.verts,3*source.nv*sizeof(float)) ||
            memcmp(prepared.faces,source.faces,3*source.nf*sizeof(int32_t))) goto done;
        for (size_t k = 0; k < 2*source.nv; k++) uv[k] = MeshBin_uv(&prepared,k);
        MeshBin_dispose(&prepared);
        fprintf(stderr,"[ink_closeup] retained-geometry UV seed verified: %s\n",prepared_uv);
    } else {
        ic_reparameterize_mode = 0;
        if (ic_reparameterize_uv(&source,uv,&pre)) goto done;
    }
    if (ic_fixed_append(&source,uv,side[0],side[1],side[2],&result,&fill,provisional?&proposal:NULL)) goto done;
    appended_nv = result.mesh.nv-source.nv; appended_nf = result.mesh.nf-source.nf;
    free(uv); uv = NULL;
    ic_reparameterize_mode = mode;
    if (c->fiber_constraints) {
        directions = calloc(result.mesh.nf,sizeof *directions);
        if (!directions || ic_measure_directions(c,&result.mesh,result.nrm,result.frame,directions,&measured)) goto done;
    }
    if (mode == 2 && ic_fixed_cut(&result,directions,&cones,&cut_edges,&fiber_cones)) goto done;
    uv = malloc(2*result.mesh.nv*sizeof(double)); if (!uv) goto done;
    if (ic_reparameterize_uv_directed(&result.mesh,directions,uv,&final)) goto done;
    (void)ic_material_measure(&result.mesh,uv,result.material,&contact);
    if (provisional) {
        patch_pass = calloc(result.npatches ? result.npatches : 1,1);
        if (!patch_pass || CloseupSurface_audit_proposal_metric(result.mesh.verts,result.mesh.nv,uv,result.frame,
                result.mesh.faces,result.mesh.nf,result.face_patch,result.npatches,patch_pass,&metric)) goto done;
        fprintf(stderr,"[ink_closeup] final proposal metric: %zu/%zu patches pass; %zu metric / %zu orientation failures\n",
                metric.passed_patches,metric.patches,metric.metric_failures,metric.orientation_failures);
    }
    snprintf(target,sizeof target,"%s/%s",c->out_dir,c->id);
    snprintf(path,sizeof path,"%s_surface.vmesh",target); ves_ensure_parent_dir(path);
    if (MeshBin_write_uv64(path,result.mesh.verts,result.mesh.nv,result.mesh.faces,result.mesh.nf,uv)) goto done;
    if (ic_material_write_aliases(target,result.mesh.nv,result.material)) goto done;
    void *output[3] = {result.nrm,result.attr,result.chart};
    for (int k = 0; k < 3; k++) {
        snprintf(path,sizeof path,"%s%s",target,suffix[k]);
        if (ic_write_file(path,output[k],result.mesh.nv*stride[k])) goto done;
    }
    snprintf(path,sizeof path,"%s_surface_source_vertex.i32",target);
    if (ic_write_file(path,result.source,result.mesh.nv*sizeof(int32_t))) goto done;
    snprintf(path,sizeof path,"%s_surface_frame_uv.f64",target);
    if (ic_write_file(path,result.frame,2*result.mesh.nv*sizeof(double))) goto done;
    if (directions) {
        snprintf(path,sizeof path,"%s_surface_fiber_arms.f64",target);
        if (ic_write_file(path,directions,result.mesh.nf*sizeof *directions)) goto done;
        if (ic_freeze_fiber_identity(c,target)) goto done;
    }
    if (provisional) {
        snprintf(path,sizeof path,"%s_surface_face_patch.i32",target);
        if (ic_write_file(path,result.face_patch,result.mesh.nf*sizeof(int32_t))) goto done;
        snprintf(path,sizeof path,"%s_surface_patch_metric_pass.u8",target);
        if (ic_write_file(path,patch_pass,result.npatches)) goto done;
    }
    snprintf(path,sizeof path,"%s_surface_report.json",target);
    report = fopen(path,"wb"); if (!report) goto done;
    ic_json_str(json_source,sizeof json_source,c->surface_source);
    fprintf(report,"{\n");
    ic_material_report(report,&contact);
    ic_growth_report(report,&final);
    if (directions) ic_fiber_report(report,&final,&measured,c);
    fprintf(report," \"schema\": \"ink-closeup-reparameterize-v1\",\n \"method\": \"slim_full_mesh\",\n"
            " \"surface_source\": %s, \"world_geometry_preserved\": true,\n"
            " \"source_vertices\": %zu, \"source_faces\": %zu,\n"
            " \"appended_vertices\": %zu, \"appended_faces\": %zu,\n"
            " \"cone_terminals\": %zu, \"fiber_cone_terminals\": %zu, \"cut_edges\": %zu,\n"
            " \"surface\": {\"verts\": %zu, \"faces\": %zu},\n"
            " \"completion\": {\"holes\": %zu, \"filled\": %zu, \"refused\": %zu},\n"
            " \"completion_provisional\": %s, \"proposal_metric_qualified\": %s,\n"
            " \"proposal_contact_refused\": %zu,\n"
            " \"proposal_metric\": {\"patches\": %zu, \"passed\": %zu, \"refused\": %zu, \"triangles\": %zu,"
            " \"metric_failures\": %zu, \"orientation_failures\": %zu, \"invalid_faces\": %zu},\n"
            " \"prepared_uv_reused\": %s, \"pre_energy_before\": %.12g, \"pre_energy_after\": %.12g,\n"
            " \"iterations\": %zu, \"accepted_steps\": %zu,\n"
            " \"energy_before\": %.12g, \"energy_after\": %.12g,\n"
            " \"flips_after\": %zu, \"boundary_intersections_after\": %zu, \"seconds\": %.3f\n}\n",
            json_source,source.nv,source.nf,appended_nv,appended_nf,cones,fiber_cones,cut_edges,result.mesh.nv,result.mesh.nf,
            fill.variational_holes,fill.variational_filled,fill.variational_refused,provisional?"true":"false",
            provisional ? (metric.metric_qualified?"true":"false") : "null",
            proposal.contact_refused,
            metric.patches,metric.passed_patches,metric.refused_patches,metric.triangles,
            metric.metric_failures,metric.orientation_failures,metric.invalid_faces,prepared_uv?"true":"false",
            pre.energy_before,pre.energy_after,final.iterations_run,final.accepted_steps,
            final.energy_before,final.energy_after,final.flips_after,final.boundary_intersections_after,
            pre.total_seconds+final.total_seconds);
    rc = fclose(report) != 0 || !ic_material_qualified(&contact) || (provisional && !metric.metric_qualified) ||
         (directions && !ic_fiber_qualified(&final)); report = NULL;
done:
    ic_reparameterize_mode = mode;
    if (report) fclose(report);
    if (rc) fprintf(stderr,"ERROR: completed-mesh repair failed; original material retained\n");
    for (int k = 0; k < 3; k++) free(side[k]);
    free(uv); free(patch_pass); free(directions); MeshBin_dispose(&source); MeshBin_dispose(&prepared); ic_fixed_free(&result);
    return rc;
}
