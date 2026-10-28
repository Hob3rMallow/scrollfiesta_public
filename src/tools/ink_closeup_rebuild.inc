/* Distinct repair contract: observed/adopted source vertices and fully observed
 * source faces are immutable mapped prefixes. Earlier generated material is
 * replaced only in this new draft; its complete source artifact stays on disk. */
typedef struct IcRebuildCounts {
    size_t observed_vertices, observed_faces, replaced_fill_vertices, replaced_inferred_faces;
    size_t retained_orphan_observed, replaced_orphan_fill;
    CloseupRebuildDomainStats domain;
} IcRebuildCounts;

static int ic_rebuild_material(const MeshBinData *source, const float *normal,
                                const uint8_t *attr, const int32_t *chart,
                                IcFixedCompletion *result, int32_t **face_source,
                                IcRebuildCounts *counts, CloseupRebuildStats *stats)
{
    CloseupGrid grid = {0}; CloseupRebuildResult proposal = {0};
    int32_t *parent = NULL, *node = NULL, *vertex = NULL, *source_to_new = NULL;
    int32_t *observed_faces = NULL, *queue = NULL;
    uint8_t *active = NULL, *domain = NULL, *outside = NULL;
    double lo[2] = {INFINITY,INFINITY}, hi[2] = {-INFINITY,-INFINITY};
    size_t nn = 0, observed_nf = 0;
    const char *phase = "source obligations";
    int rc = -1;
    memset(counts,0,sizeof *counts); *face_source = NULL;
    if (!source->nv || !source->nf || source->nv > INT32_MAX || source->nf > INT32_MAX || (!source->uv && !source->uv64)) goto done;
    parent = malloc(source->nv*sizeof *parent); node = malloc(source->nv*sizeof *node);
    source_to_new = malloc(source->nv*sizeof *source_to_new); active = calloc(source->nv,1);
    observed_faces = malloc(3*source->nf*sizeof *observed_faces);
    if (!parent || !node || !source_to_new || !active || !observed_faces) goto done;
    for (size_t v = 0; v < source->nv; v++) {
        parent[v] = (int32_t)v; node[v] = -1; source_to_new[v] = -1;
        if (attr[v] != 1 && attr[v] != 2 && attr[v] != 3) goto done;
        if (attr[v] != 2) source_to_new[v] = (int32_t)counts->observed_vertices++;
        else counts->replaced_fill_vertices++;
        for (int d = 0; d < 3; d++) if (!isfinite(source->verts[3*v+(size_t)d]) || !isfinite(normal[3*v+(size_t)d])) goto done;
        for (int d = 0; d < 2; d++) if (!isfinite(MeshBin_uv(source,2*v+(size_t)d))) goto done;
    }
    for (size_t f = 0; f < source->nf; f++) {
        int observed = 1;
        for (int k = 0; k < 3; k++) {
            int32_t v = source->faces[3*f+(size_t)k];
            if (v < 0 || (size_t)v >= source->nv) goto done;
            active[v] = 1; observed &= attr[v] != 2;
        }
        counts->observed_faces += observed != 0; counts->replaced_inferred_faces += !observed;
        for (int k = 1; k < 3; k++) {
            int32_t a = ic_fixed_root(parent,source->faces[3*f]), b = ic_fixed_root(parent,source->faces[3*f+(size_t)k]);
            if (a != b) parent[b] = a;
        }
    }
    for (size_t v = 0; v < source->nv; v++) {
        if (!active[v]) {
            counts->retained_orphan_observed += attr[v] != 2;
            counts->replaced_orphan_fill += attr[v] == 2; continue;
        }
        parent[v] = ic_fixed_root(parent,(int32_t)v);
        for (size_t k = 0; k < 2; k++) {
            double value = MeshBin_uv(source,2*v+k); lo[k] = fmin(lo[k],value); hi[k] = fmax(hi[k],value);
        }
    }
    phase = "bounded incident-material footprint";
    if (CloseupGrid_alloc(&grid,0,hi[0]-lo[0]+2*CLOSEUP_GRID_PX,
                              0,hi[1]-lo[1]+2*CLOSEUP_GRID_PX,CLOSEUP_GRID_PX)) goto done;
    grid.u0 = lo[0]-grid.h; grid.v0 = lo[1]-grid.h;
    if (grid.ni > INT32_MAX/grid.nj) goto done;
    nn = grid.ni*grid.nj;
    vertex = malloc(nn*sizeof *vertex); queue = malloc(nn*sizeof *queue);
    domain = calloc(nn,1); outside = calloc(nn,1);
    if (!vertex || !queue || !domain || !outside) goto done;
    for (size_t k = 0; k < nn; k++) vertex[k] = -1;
    for (size_t v = 0; v < source->nv; v++) if (active[v]) {
        double x = (MeshBin_uv(source,2*v)-grid.u0)/grid.h, y = (MeshBin_uv(source,2*v+1)-grid.v0)/grid.h;
        double ix = floor(x+.5), iy = floor(y+.5);
        if (fabs(x-ix)>1e-5 || fabs(y-iy)>1e-5 || ix<0 || iy<0 || ix>=(double)grid.ni || iy>=(double)grid.nj) goto done;
        size_t at = (size_t)iy*grid.ni+(size_t)ix;
        if (domain[at]) goto done; /* Ambiguous incident material cannot be collapsed into a grid sample. */
        node[v] = (int32_t)at; domain[at] = 1;
        if (attr[v] != 2) vertex[at] = source_to_new[v];
        memcpy(grid.pos+3*at,source->verts+3*v,3*sizeof(float)); memcpy(grid.nrm+3*at,normal+3*v,3*sizeof(float));
        grid.state[at] = attr[v] == 2 ? CLOSEUP_NODE_FILL : CLOSEUP_NODE_DATA;
        grid.chart[at] = chart[v]; grid.component[at] = parent[v];
    }
    /* One padded exterior flood, with the same eight-neighbor convention as
     * native enclosed-gap detection. No apron is added outside this footprint. */
    size_t tail = 1;
    queue[0] = 0; outside[0] = 1;
    for (size_t head = 0; head < tail; head++) {
        size_t at = (size_t)queue[head], x = at%grid.ni, y = at/grid.ni;
        for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
            if ((!dx && !dy) || (dx<0 && !x) || (dy<0 && !y) ||
                (dx>0 && x+1>=grid.ni) || (dy>0 && y+1>=grid.nj)) continue;
            size_t next = (size_t)((int64_t)at+(int64_t)dy*(int64_t)grid.ni+dx);
            if (domain[next] || outside[next]) continue;
            outside[next] = 1; queue[tail++] = (int32_t)next;
        }
    }
    for (size_t k = 0; k < nn; k++) if (!outside[k]) domain[k] = 1;
    for (size_t f = 0; f < source->nf; f++) {
        int observed = 1;
        for (int k = 0; k < 3; k++) observed &= attr[source->faces[3*f+(size_t)k]] != 2;
        if (!observed) continue;
        for (int k = 0; k < 3; k++) observed_faces[3*observed_nf+(size_t)k] = node[source->faces[3*f+(size_t)k]];
        observed_nf++;
    }
    if (observed_nf != counts->observed_faces) goto done;
    phase = "bounded DATA-supported regrowth proposal domain";
    if (CloseupSurface_rebuild_regrowth_domain(&grid,domain,CLOSEUP_COMPLETE_FILL_PX,&counts->domain)) goto done;
    fprintf(stderr,"[ink_closeup] supported-regrowth-domain-v2: %zu required + %zu added = %zu proposed nodes; "
            "%zu incident DATA donors, %.9g frame-pixel reach, %zu required beyond reach\n",
            counts->domain.required_nodes,counts->domain.added_nodes,counts->domain.nodes,counts->domain.data_donors,
            counts->domain.reach_px,counts->domain.required_beyond_reach);
    phase = "DATA-anchored generated-material proposal";
    if (CloseupSurface_rebuild_generated(&grid,domain,observed_faces,observed_nf,CLOSEUP_COMPLETE_FILL_PX,&proposal,stats)) goto done;
    if (counts->observed_vertices+proposal.nodes > INT32_MAX || counts->observed_faces+proposal.faces > SIZE_MAX/(3*sizeof(int32_t))) goto done;
    result->mesh.nv = counts->observed_vertices+proposal.nodes;
    result->mesh.nf = counts->observed_faces+proposal.faces; result->npatches = proposal.patches;
    result->mesh.verts = malloc(3*result->mesh.nv*sizeof(float)); result->mesh.uv64 = malloc(2*result->mesh.nv*sizeof(double));
    result->mesh.faces = malloc(3*result->mesh.nf*sizeof(int32_t)); result->nrm = malloc(3*result->mesh.nv*sizeof(float));
    result->attr = malloc(result->mesh.nv); result->chart = malloc(result->mesh.nv*sizeof(int32_t));
    result->source = malloc(result->mesh.nv*sizeof(int32_t)); result->frame = malloc(2*result->mesh.nv*sizeof(double));
    result->face_patch = malloc(result->mesh.nf*sizeof(int32_t)); *face_source = malloc(result->mesh.nf*sizeof(int32_t));
    if (!result->mesh.verts || !result->mesh.uv64 || !result->mesh.faces || !result->nrm || !result->attr ||
        !result->chart || !result->source || !result->frame || !result->face_patch || !*face_source) goto done;
    phase = "exact observed prefixes";
    for (size_t v = 0; v < source->nv; v++) if (source_to_new[v] >= 0) {
        size_t to = (size_t)source_to_new[v];
        memcpy(result->mesh.verts+3*to,source->verts+3*v,3*sizeof(float)); memcpy(result->nrm+3*to,normal+3*v,3*sizeof(float));
        result->attr[to] = attr[v]; result->chart[to] = chart[v]; result->source[to] = (int32_t)v;
        for (size_t d = 0; d < 2; d++) result->frame[2*to+d] = result->mesh.uv64[2*to+d] = MeshBin_uv(source,2*v+d);
    }
    size_t next_vertex = counts->observed_vertices, next_face = 0;
    for (size_t k = 0; k < nn; k++) if (proposal.node_patch[k] >= 0) {
        size_t to = next_vertex++; vertex[k] = (int32_t)to;
        memcpy(result->mesh.verts+3*to,grid.pos+3*k,3*sizeof(float)); memcpy(result->nrm+3*to,grid.nrm+3*k,3*sizeof(float));
        result->attr[to] = 2; result->chart[to] = -1; result->source[to] = -1;
        result->frame[2*to] = result->mesh.uv64[2*to] = grid.u0+grid.h*(double)(k%grid.ni);
        result->frame[2*to+1] = result->mesh.uv64[2*to+1] = grid.v0+grid.h*(double)(k/grid.ni);
    }
    if (next_vertex != result->mesh.nv) goto done;
    for (size_t f = 0; f < source->nf; f++) {
        int observed = 1;
        for (int k = 0; k < 3; k++) observed &= source_to_new[source->faces[3*f+(size_t)k]] >= 0;
        if (!observed) continue;
        for (int k = 0; k < 3; k++) result->mesh.faces[3*next_face+(size_t)k] = source_to_new[source->faces[3*f+(size_t)k]];
        (*face_source)[next_face] = (int32_t)f; result->face_patch[next_face++] = -1;
    }
    if (next_face != counts->observed_faces) goto done;
    for (size_t f = 0; f < proposal.faces; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t to = vertex[proposal.face_nodes[3*f+(size_t)k]];
            if (to < 0) goto done;
            result->mesh.faces[3*next_face+(size_t)k] = to;
        }
        (*face_source)[next_face] = -1; result->face_patch[next_face++] = proposal.face_patch[f];
    }
    fprintf(stderr,"[ink_closeup] DATA rebuild prefixes: %zu observed vertices (%zu orphans), %zu observed faces; "
            "replaced %zu FILL vertices (%zu orphans), %zu inferred faces; original artifact retained\n",
            counts->observed_vertices,counts->retained_orphan_observed,counts->observed_faces,
            counts->replaced_fill_vertices,counts->replaced_orphan_fill,counts->replaced_inferred_faces);
    rc = 0;
done:
    if (rc) { fprintf(stderr,"ERROR: DATA-preserving rebuild failed at %s\n",phase); ic_fixed_free(result); free(*face_source); *face_source = NULL; }
    CloseupGrid_free(&grid); CloseupRebuildResult_free(&proposal);
    free(parent); free(node); free(vertex); free(source_to_new); free(observed_faces); free(queue); free(active); free(domain); free(outside);
    return rc;
}

static int ic_rebuild_fresh_target(const char *target)
{
    const char *suffix[] = {"_surface.vmesh","_surface_nrm.f32","_surface_attr.u8","_surface_chart.i32",
        "_surface_source_vertex.i32","_surface_source_face.i32","_surface_frame_uv.f64",
        "_surface_face_patch.i32","_surface_patch_metric_pass.u8","_surface_fiber_arms.f64",
        "_surface_inference_identity.json","_surface_material_vertex.i32","_surface_report.json"};
    char path[2800];
    for (size_t k = 0; k < sizeof suffix/sizeof suffix[0]; k++) {
        snprintf(path,sizeof path,"%s%s",target,suffix[k]);
        FILE *existing = fopen(path,"rb");
        if (existing) { fclose(existing); fprintf(stderr,"ERROR: rebuild output already exists: %s\n",path); return -1; }
    }
    return 0;
}

static int ic_rebuild_world_qualified(const CloseupRebuildStats *st)
{
    return st->world_exact.complete && st->world_proximity.complete &&
           !st->world_exact.observed_observed && !st->world_proximity.observed_observed;
}

static void ic_rebuild_world_report(FILE *report, const CloseupRebuildStats *st)
{
    fprintf(report," \"world_geometry_contacts\": {\"policy\": \"IntersectionCleanup_complete_readonly_v1\","
        "\"scope\": \"all_staged_faces_before_atomic_patch_rejection\",\"all_faces_queried\": true,"
        "\"inferred_conflicts_rejected_as_whole_patches\": true,\"original_observed_faces_retained\": true,"
        "\"retained_inferred_conflicts\": 0,\"world_geometry_qualified\": %s,",
        ic_rebuild_world_qualified(st) ? "true" : "false");
    for (int mode = 0; mode < 2; mode++) {
        const CloseupRebuildWorldStats *s = mode ? &st->world_proximity : &st->world_exact;
        fprintf(report,"\"%s\":{\"gap_vox\":%d,\"parallel_angle_deg\":20,\"complete\":%s,\"pairs\":%zu,"
            "\"retained_observed_retained_observed\":%zu,\"retained_observed_inferred\":%zu,\"inferred_inferred\":%zu,"
            "\"overlaps\":%zu,\"stabs\":%zu,\"folds\":%zu}%s",
            mode ? "near_coplanar" : "exact",mode,s->complete ? "true" : "false",s->pairs,s->observed_observed,
            s->observed_inferred,s->inferred_inferred,s->overlaps,s->stabs,s->folds,mode ? "},\n" : ",");
    }
}

static int ic_stage_rebuild_completed(const IcConfig *c)
{
    MeshBinData source = {0}; IcFixedCompletion result = {0}; IcRebuildCounts counts = {0};
    CloseupRebuildStats rebuild = {0}; CloseupProposalMetricStats metric = {0};
    SlimRefineStats final = {0}; CloseupFiberStats measured = {0}; SlimRefineDirection *directions = NULL;
    AsmContactStats contact = {0};
    int32_t *face_source = NULL; uint8_t *patch_pass = NULL; double *uv = NULL;
    const char *suffix[3] = {"_surface_nrm.f32","_surface_attr.u8","_surface_chart.i32"};
    size_t stride[3] = {3*sizeof(float),1,sizeof(int32_t)}, got = 0, cones = 0, cut_edges = 0, fiber_cones = 0;
    void *side[3] = {0}; char path[2800], target[2600], json_source[4800];
    FILE *report = NULL; int rc = 1;
    snprintf(target,sizeof target,"%s/%s",c->out_dir,c->id);
    if (ic_rebuild_fresh_target(target)) goto done;
    snprintf(path,sizeof path,"%s_surface.vmesh",c->surface_source);
    if (MeshBin_read_precise_malloc(path,&source)) goto done;
    for (int k = 0; k < 3; k++) {
        snprintf(path,sizeof path,"%s%s",c->surface_source,suffix[k]);
        side[k] = ic_read_file(path,source.nv*stride[k],&got);
        if (!side[k] || got != source.nv*stride[k]) goto done;
    }
    if (ic_rebuild_material(&source,side[0],side[1],side[2],&result,&face_source,&counts,&rebuild)) goto done;
    /* Measure on real material neighborhoods before vertex-only cone cuts.
     * Faces keep their order, so the physical directions remain correctly bound. */
    if (c->fiber_constraints) {
        directions = calloc(result.mesh.nf,sizeof *directions);
        if (!directions || ic_measure_directions(c,&result.mesh,result.nrm,result.frame,directions,&measured)) goto done;
    }
    if (ic_reparameterize_mode == 2 && ic_fixed_cut(&result,directions,&cones,&cut_edges,&fiber_cones)) goto done;
    uv = malloc(2*result.mesh.nv*sizeof(double)); patch_pass = calloc(result.npatches ? result.npatches : 1,1);
    if (!uv || !patch_pass || ic_reparameterize_uv_directed(&result.mesh,directions,uv,&final)) goto done;
    (void)ic_material_measure(&result.mesh,uv,result.material,&contact);
    if (CloseupSurface_audit_proposal_metric(result.mesh.verts,result.mesh.nv,uv,result.frame,result.mesh.faces,result.mesh.nf,
            result.face_patch,result.npatches,patch_pass,&metric)) goto done;
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
    snprintf(path,sizeof path,"%s_surface_source_face.i32",target);
    if (ic_write_file(path,face_source,result.mesh.nf*sizeof(int32_t))) goto done;
    snprintf(path,sizeof path,"%s_surface_frame_uv.f64",target);
    if (ic_write_file(path,result.frame,2*result.mesh.nv*sizeof(double))) goto done;
    snprintf(path,sizeof path,"%s_surface_face_patch.i32",target);
    if (ic_write_file(path,result.face_patch,result.mesh.nf*sizeof(int32_t))) goto done;
    snprintf(path,sizeof path,"%s_surface_patch_metric_pass.u8",target);
    if (ic_write_file(path,patch_pass,result.npatches)) goto done;
    if (directions) {
        snprintf(path,sizeof path,"%s_surface_fiber_arms.f64",target);
        if (ic_write_file(path,directions,result.mesh.nf*sizeof *directions)) goto done;
        if (ic_freeze_fiber_identity(c,target)) goto done;
    }
    snprintf(path,sizeof path,"%s_surface_report.json",target); report = fopen(path,"wb"); if (!report) goto done;
    ic_json_str(json_source,sizeof json_source,c->surface_source);
    fprintf(report,"{\n");
    ic_material_report(report,&contact);
    ic_rebuild_world_report(report,&rebuild);
    ic_growth_report(report,&final);
    if (directions) ic_fiber_report(report,&final,&measured,c);
    fprintf(report," \"schema\": \"ink-closeup-data-rebuild-v1\",\n \"method\": \"slim_full_mesh\",\n"
        " \"preservation_schema\": \"closeup-data-preserving-rebuild-v1\",\n"
        " \"surface_source\": %s, \"source_artifact_modified\": false, \"observed_geometry_preserved\": true,\n"
        " \"world_geometry_preserved\": true, \"world_geometry_reference\": \"staged_data_preserving_rebuild\",\n"
        " \"source_vertices\": %zu, \"source_faces\": %zu,\n"
        " \"observed_vertex_prefix\": %zu, \"observed_face_prefix\": %zu,\n"
        " \"retained_orphan_observed_vertices\": %zu, \"replaced_generated_vertices\": %zu,"
        " \"replaced_orphan_generated_vertices\": %zu, \"replaced_inferred_faces\": %zu,\n"
        " \"generated_vertices\": %zu, \"generated_faces\": %zu, \"cone_terminals\": %zu, \"fiber_cone_terminals\": %zu, \"cut_edges\": %zu,\n"
        " \"surface\": {\"verts\": %zu, \"faces\": %zu},\n"
        " \"completion_provisional\": true, \"proposal_metric_qualified\": %s,\n"
        " \"completion\": {\"patches\": %zu, \"proposed\": %zu, \"refused\": %zu,"
        " \"reach_refused\": %zu, \"orientation_refused\": %zu, \"contact_refused\": %zu,"
        " \"unsupported_nodes\": %zu, \"separated_cells\": %zu, \"alternate_inferred_diagonals\": %zu},\n"
        " \"domain\": {\"policy\": \"supported-regrowth-domain-v2\", \"required_policy\": \"prior_incident_material_plus_enclosed_gaps\","
        " \"required_nodes\": %zu, \"required_enclosed_gap_nodes\": %zu, \"added_nodes\": %zu, \"nodes\": %zu,"
        " \"required_beyond_reach\": %zu, \"data_donors\": %zu, \"reach_px\": %.12g,"
        " \"source_frame_bounds\": [%.12g,%.12g,%.12g,%.12g], \"optional_added_nodes\": true},\n"
        " \"component_support\": {\"source_components\": %zu, \"supported_groups\": %zu, \"certified_joins\": %zu, \"refused_joins\": %zu},\n"
        " \"geometry_screen\": {\"policy\": \"immutable_DATA_tangent\", \"weight\": \"(h2/(h2+frame_distance_squared))^2\","
        " \"screened_attempt_nodes\": %zu, \"screen_energy\": %.12g},\n"
        " \"orientation_certificate\": {\"policy\": \"observed_anchored_oriented_connectivity_and_world_contact_v2\","
        " \"normal_angle_is_diagnostic\": true,"
        " \"components\": %zu, \"observed_anchored_components\": %zu, \"unproven_components\": %zu, \"unproven_faces\": %zu,"
        " \"degenerate_faces\": %zu, \"nonpositive_normal_edges\": %zu, \"winding_conflicts\": %zu, \"nonmanifold_edges\": %zu},\n"
        " \"generated_motion\": {\"initialization_policy\": \"DATA_variational_target_no_distant_normal_projection\","
        " \"old_fill_is_initialization_only\": true, \"old_fill_start_nodes\": %zu,"
        " \"target_start_patches\": %zu, \"sweeps\": %zu, \"vertex_moves\": %zu, \"invalid_before\": %zu, \"invalid_after\": %zu,"
        " \"target_error_before\": %.12g, \"target_error_after\": %.12g,"
        " \"coupled_passes\": %zu, \"coupled_blocks\": %zu, \"coupled_vertex_moves\": %zu, \"coupled_TAUCS_solves\": %zu},\n"
        " \"proposal_metric\": {\"patches\": %zu, \"passed\": %zu, \"refused\": %zu, \"triangles\": %zu,"
        " \"metric_failures\": %zu, \"orientation_failures\": %zu, \"invalid_faces\": %zu},\n"
        " \"iterations\": %zu, \"accepted_steps\": %zu, \"energy_before\": %.12g, \"energy_after\": %.12g,\n"
        " \"flips_after\": %zu, \"boundary_intersections_after\": %zu, \"seconds\": %.3f, \"final_qualified\": false\n}\n",
        json_source,source.nv,source.nf,counts.observed_vertices,counts.observed_faces,counts.retained_orphan_observed,
        counts.replaced_fill_vertices,counts.replaced_orphan_fill,counts.replaced_inferred_faces,rebuild.proposed_nodes,rebuild.proposed_faces,
        cones,fiber_cones,cut_edges,result.mesh.nv,result.mesh.nf,metric.metric_qualified ? "true" : "false",
        rebuild.patches,rebuild.proposed_patches,rebuild.refused_patches,rebuild.reach_refused,rebuild.orientation_refused,rebuild.contact_refused,
        rebuild.unsupported_nodes,rebuild.separated_cells,rebuild.alternate_diagonals,
        counts.domain.required_nodes,counts.domain.required_enclosed_gap_nodes,counts.domain.added_nodes,counts.domain.nodes,
        counts.domain.required_beyond_reach,counts.domain.data_donors,counts.domain.reach_px,
        counts.domain.source_frame_bounds[0],counts.domain.source_frame_bounds[1],counts.domain.source_frame_bounds[2],counts.domain.source_frame_bounds[3],
        rebuild.source_components,rebuild.supported_groups,rebuild.certified_joins,rebuild.refused_joins,
        rebuild.screened_nodes,rebuild.tangent_screen_energy,
        rebuild.transport_components,rebuild.transport_anchored_components,rebuild.transport_unanchored_components,rebuild.transport_unanchored_faces,
        rebuild.transport_degenerate_faces,rebuild.transport_bad_edges,rebuild.transport_winding_edges,rebuild.transport_nonmanifold_edges,
        rebuild.old_fill_start_nodes,rebuild.target_start_patches,rebuild.motion_sweeps,rebuild.motion_vertices,rebuild.orientation_before_motion,
        rebuild.orientation_after_motion,rebuild.target_error_before_motion,rebuild.target_error_after_motion,
        rebuild.coupled_passes,rebuild.coupled_blocks,rebuild.coupled_vertices,rebuild.coupled_solves,
        metric.patches,metric.passed_patches,metric.refused_patches,metric.triangles,metric.metric_failures,metric.orientation_failures,metric.invalid_faces,
        final.iterations_run,final.accepted_steps,final.energy_before,final.energy_after,final.flips_after,final.boundary_intersections_after,final.total_seconds);
    rc = fclose(report) != 0 || !ic_material_qualified(&contact) || !ic_rebuild_world_qualified(&rebuild) || !metric.metric_qualified ||
         (directions && !ic_fiber_qualified(&final)); report = NULL;
done:
    if (report) fclose(report);
    if (rc) fprintf(stderr,"ERROR: DATA rebuild did not qualify; original artifact retained, any new output is diagnostic only\n");
    for (int k = 0; k < 3; k++) free(side[k]);
    free(uv); free(patch_pass); free(face_source); free(directions); MeshBin_dispose(&source); ic_fixed_free(&result);
    return rc;
}

static int ic_selftest_rebuild_material(void)
{
    int fails = 0;
    for (int reverse = 0; reverse < 2; reverse++) {
        CloseupGrid grid = {0}; MeshBinData source = {0}; IcFixedCompletion result = {0};
        CloseupRebuildStats stats = {0}; IcRebuildCounts counts = {0};
        int32_t *nodes = NULL, *chart = NULL, *face_source = NULL;
        float *normal = NULL; uint8_t *attr = NULL; double *grid_uv = NULL;
        int before = fails;
        if (CloseupGrid_alloc(&grid,-12,12,-12,12,2)) return 1;
        size_t nn = grid.ni*grid.nj;
        grid_uv = malloc(2*nn*sizeof(double));
        if (!grid_uv) { CloseupGrid_free(&grid); return 1; }
        for (size_t k = 0; k < nn; k++) {
            double u = grid.u0+grid.h*(double)(k%grid.ni), v = grid.v0+grid.h*(double)(k/grid.ni);
            grid.state[k] = fabs(u)<6 && fabs(v)<6 ? CLOSEUP_NODE_FILL : CLOSEUP_NODE_DATA;
            grid.component[k] = 7; grid.chart[k] = grid.state[k] == CLOSEUP_NODE_DATA ? 19 : -1;
            grid.pos[3*k] = (float)(12000+u); grid.pos[3*k+1] = (float)(8000+v);
            grid.pos[3*k+2] = grid.state[k] == CLOSEUP_NODE_DATA ? 4000 : 4030;
            grid.nrm[3*k+2] = reverse ? -1 : 1;
            grid_uv[2*k] = u; grid_uv[2*k+1] = v;
        }
        if (CloseupSurface_mesh(&grid,grid_uv,&source.verts,&source.uv64,&normal,&source.faces,&nodes,&source.nv,&source.nf)) {
            fails++; goto test_done;
        }
        size_t original_nv = source.nv;
        float *new_xyz = realloc(source.verts,3*(original_nv+2)*sizeof(float));
        if (!new_xyz) { fails++; goto test_done; } source.verts = new_xyz;
        double *new_uv = realloc(source.uv64,2*(original_nv+2)*sizeof(double));
        if (!new_uv) { fails++; goto test_done; } source.uv64 = new_uv;
        float *new_normal = realloc(normal,3*(original_nv+2)*sizeof(float));
        if (!new_normal) { fails++; goto test_done; } normal = new_normal;
        source.nv += 2;
        attr = malloc(source.nv); chart = malloc(source.nv*sizeof(int32_t));
        if (!attr || !chart) { fails++; goto test_done; }
        for (size_t v = 0; v < original_nv; v++) { attr[v] = grid.state[nodes[v]]; chart[v] = grid.chart[nodes[v]]; }
        attr[0] = 3; /* adopted observations retain their provenance value */
        for (size_t k = original_nv; k < source.nv; k++) {
            source.verts[3*k] = 12000; source.verts[3*k+1] = 8000; source.verts[3*k+2] = 9000;
            source.uv64[2*k] = k == original_nv ? 0 : 200; source.uv64[2*k+1] = 0;
            normal[3*k] = normal[3*k+1] = 0; normal[3*k+2] = reverse ? -1 : 1;
            attr[k] = k == original_nv ? 3 : 2; chart[k] = k == original_nv ? 23 : -1;
        }
        if (reverse) for (size_t f = 0; f < source.nf; f++) {
            int32_t swap = source.faces[3*f+1]; source.faces[3*f+1] = source.faces[3*f+2]; source.faces[3*f+2] = swap;
        }
        if (ic_rebuild_material(&source,normal,attr,chart,&result,&face_source,&counts,&stats)) { fails++; goto test_done; }
        size_t next = 0, observed_nf = 0;
        int32_t orphan_copy = -1;
        for (size_t v = 0; v < source.nv; v++) if (attr[v] != 2) {
            if (result.source[next] != (int32_t)v || result.attr[next] != attr[v] || result.chart[next] != chart[v] ||
                memcmp(result.mesh.verts+3*next,source.verts+3*v,3*sizeof(float)) ||
                memcmp(result.nrm+3*next,normal+3*v,3*sizeof(float)) ||
                memcmp(result.frame+2*next,source.uv64+2*v,2*sizeof(double))) fails++;
            if (v == original_nv) orphan_copy = (int32_t)next;
            next++;
        }
        for (size_t f = 0; f < source.nf; f++) {
            int observed = 1;
            for (int k = 0; k < 3; k++) observed &= attr[source.faces[3*f+(size_t)k]] != 2;
            if (!observed) continue;
            if (face_source[observed_nf] != (int32_t)f || result.face_patch[observed_nf] != -1) fails++;
            for (int k = 0; k < 3; k++) if (result.source[result.mesh.faces[3*observed_nf+(size_t)k]] != source.faces[3*f+(size_t)k]) fails++;
            observed_nf++;
        }
        for (size_t k = 0; k < 3*result.mesh.nf; k++) if (result.mesh.faces[k] == orphan_copy) fails++;
        if (counts.retained_orphan_observed != 1 || counts.replaced_orphan_fill != 1 ||
            counts.replaced_fill_vertices != stats.old_fill_nodes+1 || stats.refused_patches ||
            counts.observed_vertices != next || counts.observed_faces != observed_nf || !stats.proposed_nodes) fails++;
test_done:
        fprintf(stderr,"  DATA rebuild driver winding %d: observed/adopted prefixes and inactive orphans %s\n",reverse,fails == before ? "ok" : "FAIL");
        free(nodes); free(chart); free(face_source); free(normal); free(attr); free(grid_uv);
        MeshBin_dispose(&source); ic_fixed_free(&result); CloseupGrid_free(&grid);
    }
    fprintf(stderr,"[ink_closeup] DATA rebuild driver: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails != 0;
}
