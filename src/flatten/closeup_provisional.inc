/* Proposal staging is deliberately separate from accepted fixed-grid infill.
 * A proposal may need a cut before its metric is feasible. Neither the grid
 * modification nor a successful staging return qualifies its final UV. */
int CloseupSurface_propose_fixed(CloseupGrid *g, double max_distance_px, int32_t **patch_ids,
                                 CloseupProposalStats *st)
{
    CloseupFitStats fit = {0}; CloseupDedupStats contact = {0};
    CloseupGrid audit = {0};
    int32_t *patch = NULL, *remap = NULL;
    uint8_t *refused = NULL;
    size_t nn = 0, npatches = 0;
    int rc = -1;
    if (!g || !patch_ids || !st) return -1;
    *patch_ids = NULL; memset(st,0,sizeof *st); st->metric_qualification_pending = 1;
    if (!g->ni || !g->nj || g->ni > INT32_MAX/g->nj) return -1;
    nn = g->ni*g->nj;
    patch = (int32_t *)malloc(nn*sizeof(int32_t));
    if (!patch) goto done;
    for (size_t k = 0; k < nn; k++) patch[k] = -1;
    if (cs_complete_fixed(g,NULL,max_distance_px,&fit,1,patch)) goto done;
    npatches = fit.variational_filled;
    refused = (uint8_t *)calloc(npatches ? npatches : 1,1);
    remap = (int32_t *)malloc((npatches ? npatches : 1)*sizeof(int32_t));
    if (!refused || !remap) goto done;
    /* Protect retained FILL as well as DATA in an audit-only copy. The native
     * contact rule may remove audit nodes; the source and proposal collars
     * remain untouched. One offending new node rejects its entire patch. */
    audit = *g;
    audit.state = (uint8_t *)malloc(nn);
    audit.pos = (float *)malloc(3*nn*sizeof(float));
    audit.chart = (int32_t *)malloc(nn*sizeof(int32_t));
    audit.component = (int32_t *)malloc(nn*sizeof(int32_t));
    if (!audit.state || !audit.pos || !audit.chart || !audit.component) goto done;
    memcpy(audit.state,g->state,nn); memcpy(audit.pos,g->pos,3*nn*sizeof(float));
    memcpy(audit.chart,g->chart,nn*sizeof(int32_t)); memcpy(audit.component,g->component,nn*sizeof(int32_t));
    for (size_t k = 0; k < nn; k++) if (patch[k] < 0 && audit.state[k]) audit.state[k] = CLOSEUP_NODE_DATA;
    if (CloseupSurface_dedup_fill(&audit,CLOSEUP_DEDUP_VOX,CLOSEUP_DEDUP_FRAME_PX,CLOSEUP_DEDUP_COS,&contact)) goto done;
    for (size_t k = 0; k < nn; k++) if (patch[k] >= 0 && !audit.state[k]) refused[patch[k]] = 1;
    st->holes = fit.variational_holes; st->refused = fit.variational_refused;
    for (size_t p = 0; p < npatches; p++) {
        if (refused[p]) { remap[p] = -1; st->contact_refused++; st->refused++; }
        else remap[p] = (int32_t)st->proposed++;
    }
    for (size_t k = 0; k < nn; k++) {
        if (patch[k] < 0) continue;
        patch[k] = remap[patch[k]];
        if (patch[k] >= 0) { st->nodes++; continue; }
        g->state[k] = CLOSEUP_NODE_NONE; g->chart[k] = -1; g->component[k] = -1; g->support[k] = 0;
        for (int d = 0; d < 3; d++) {
            g->pos[3*k+(size_t)d] = NAN;
            g->tu[3*k+(size_t)d] = g->tv[3*k+(size_t)d] = g->nrm[3*k+(size_t)d] = 0;
        }
    }
    fprintf(stderr,"  Provisional Cartesian infill: %zu/%zu patches staged, %zu refused (%zu contact); "
            "%zu nodes; final metric qualification PENDING\n",st->proposed,st->holes,st->refused,st->contact_refused,st->nodes);
    *patch_ids = patch; patch = NULL; rc = 0;
done:
    free(patch); free(refused); free(remap);
    /* Only these four fields belong to the shallow audit copy. */
    free(audit.state); free(audit.pos); free(audit.chart); free(audit.component);
    return rc;
}

int CloseupSurface_audit_proposal_metric(const float *verts, size_t nv, const double *uv, const double *frame_uv,
                                        const int32_t *faces, size_t nf, const int32_t *face_patch,
                                        size_t npatches, uint8_t *patch_pass, CloseupProposalMetricStats *st)
{
    size_t *counts = NULL;
    if (!st) return -1;
    memset(st,0,sizeof *st); st->patches = npatches; st->min_squared_scale = INFINITY;
    if (!verts || !uv || !frame_uv || !faces || !face_patch || (npatches && !patch_pass)) return -1;
    counts = (size_t *)calloc(npatches ? npatches : 1,sizeof(size_t));
    if (!counts) return -1;
    if (npatches) memset(patch_pass,1,npatches);
    for (size_t f = 0; f < nf; f++) {
        int32_t patch = face_patch[f];
        if (patch == -1) continue;
        if (patch < 0 || (size_t)patch >= npatches) { free(counts); return -1; }
        size_t v[3]; int valid = 1;
        st->triangles++; counts[patch]++;
        for (int k = 0; k < 3; k++) {
            int32_t q = faces[3*f+(size_t)k];
            if (q < 0 || (size_t)q >= nv) { valid = 0; break; }
            v[k] = (size_t)q;
            for (int d = 0; d < 3; d++) valid &= isfinite(verts[3*v[k]+(size_t)d]) != 0;
            for (int d = 0; d < 2; d++) valid &= isfinite(uv[2*v[k]+(size_t)d]) && isfinite(frame_uv[2*v[k]+(size_t)d]);
        }
        if (!valid) { st->invalid_faces++; patch_pass[patch] = 0; continue; }
        double au = uv[2*v[1]]-uv[2*v[0]], av = uv[2*v[1]+1]-uv[2*v[0]+1];
        double bu = uv[2*v[2]]-uv[2*v[0]], bv = uv[2*v[2]+1]-uv[2*v[0]+1];
        double det = au*bv-av*bu;
        double fu = frame_uv[2*v[1]]-frame_uv[2*v[0]], fv = frame_uv[2*v[1]+1]-frame_uv[2*v[0]+1];
        double gu = frame_uv[2*v[2]]-frame_uv[2*v[0]], gv = frame_uv[2*v[2]+1]-frame_uv[2*v[0]+1];
        double reference_det = fu*gv-fv*gu;
        if (!isfinite(det) || !isfinite(reference_det) || !(fabs(det) > 0) || !(fabs(reference_det) > 0) ||
            (det > 0) != (reference_det > 0)) {
            st->orientation_failures++; patch_pass[patch] = 0; continue;
        }
        double e[3], across[3], cross[3], ee = 0, ff = 0, ef = 0, area2 = 0;
        for (int d = 0; d < 3; d++) {
            double a = (double)verts[3*v[1]+(size_t)d]-verts[3*v[0]+(size_t)d];
            double b = (double)verts[3*v[2]+(size_t)d]-verts[3*v[0]+(size_t)d];
            e[d] = (a*bv-b*av)/det; across[d] = (b*au-a*bu)/det;
            ee += e[d]*e[d]; ff += across[d]*across[d]; ef += e[d]*across[d];
        }
        for (int d = 0; d < 3; d++) {
            cross[d] = e[(d+1)%3]*across[(d+2)%3]-e[(d+2)%3]*across[(d+1)%3];
            area2 += cross[d]*cross[d];
        }
        double hi = .5*(ee+ff+hypot(ee-ff,2*ef));
        double lo = hi > 0 ? area2/hi : 0;
        if (isfinite(lo) && lo < st->min_squared_scale) st->min_squared_scale = lo;
        if (isfinite(hi) && hi > st->max_squared_scale) st->max_squared_scale = hi;
        if (!isfinite(lo) || !isfinite(hi) || lo < .25 || hi > 4) {
            st->metric_failures++; patch_pass[patch] = 0;
        }
    }
    for (size_t p = 0; p < npatches; p++) {
        if (!counts[p]) { st->empty_patches++; patch_pass[p] = 0; }
        if (patch_pass[p]) st->passed_patches++; else st->refused_patches++;
    }
    if (!isfinite(st->min_squared_scale)) st->min_squared_scale = 0;
    st->metric_qualified = st->refused_patches == 0;
    free(counts); return 0;
}

static int cs_selftest_provisional(void)
{
    int fails = 0;
    for (int test = 0; test < 2; test++) {
        CloseupGrid g = {0}; CloseupProposalStats proposal = {0};
        int32_t *patch = NULL;
        if (CloseupGrid_alloc(&g,-8,test ? 64 : 8,-8,8,2)) return 1;
        size_t nn = g.ni*g.nj, hole = 0;
        float *original = (float *)malloc(3*nn*sizeof(float));
        uint8_t *original_state = (uint8_t *)malloc(nn);
        if (!original || !original_state) { free(original); free(original_state); CloseupGrid_free(&g); return 1; }
        for (size_t k = 0; k < nn; k++) {
            double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
            if (u > 8 && (!test || u < 58 || fabs(v)>2)) continue;
            g.state[k] = CLOSEUP_NODE_DATA; g.component[k] = u>8 ? 2 : 1;
            g.pos[3*k] = (float)(12000+.05*(u>8 ? u-60 : u));
            g.pos[3*k+1] = (float)(8000+.05*v); g.pos[3*k+2] = u>8 ? 4000.2f : 4000;
            g.tu[3*k] = .05f; g.tv[3*k+1] = .05f; cs_set_normal(&g,k);
            if (fabs(u)<3 && fabs(v)<3) {
                g.state[k] = CLOSEUP_NODE_NONE; g.component[k] = -1; hole++;
                g.pos[3*k] = g.pos[3*k+1] = g.pos[3*k+2] = NAN;
            }
        }
        memcpy(original,g.pos,3*nn*sizeof(float)); memcpy(original_state,g.state,nn);
        if (CloseupSurface_propose_fixed(&g,32,&patch,&proposal)) { fails++; goto test_done; }
        if (!proposal.metric_qualification_pending) fails++;
        for (size_t k = 0; k < nn; k++) if (original_state[k] &&
            (g.state[k] != original_state[k] || memcmp(original+3*k,g.pos+3*k,3*sizeof(float)))) fails++;
        if (test) {
            if (proposal.proposed || proposal.nodes || proposal.contact_refused != 1) fails++;
        } else {
            float *verts = NULL, *normals = NULL; double *frame = NULL, *uv = NULL;
            int32_t *faces = NULL, *nodes = NULL, *face_patch = NULL; size_t nv = 0, nf = 0;
            double *grid_uv = (double *)malloc(2*nn*sizeof(double));
            if (proposal.proposed != 1 || proposal.nodes != hole || !grid_uv) { free(grid_uv); fails++; goto test_done; }
            for (size_t k = 0; k < nn; k++) {
                grid_uv[2*k] = g.u0+g.h*(double)(k%g.ni);
                grid_uv[2*k+1] = g.v0+g.h*(double)(k/g.ni);
            }
            if (CloseupSurface_mesh(&g,grid_uv,&verts,&frame,&normals,&faces,&nodes,&nv,&nf)) {
                free(grid_uv); fails++; goto test_done;
            }
            uv = (double *)malloc(2*nv*sizeof(double)); face_patch = (int32_t *)malloc(nf*sizeof(int32_t));
            if (!uv || !face_patch) fails++;
            else {
                for (size_t f = 0; f < nf; f++) {
                    face_patch[f] = -1;
                    for (int k = 0; k < 3; k++) if (patch[nodes[faces[3*f+(size_t)k]]] >= 0) face_patch[f] = 0;
                }
                for (int trial = 0; trial < 4; trial++) {
                    uint8_t passed = 0; CloseupProposalMetricStats metric = {0};
                    for (size_t v = 0; v < nv; v++) {
                        uv[2*v] = frame[2*v]*(trial == 0 ? 1 : trial == 2 ? -.05 : .05);
                        uv[2*v+1] = frame[2*v+1]*(trial == 0 ? 1 : .05);
                    }
                    if (trial == 3) for (size_t f = 0; f < nf; f++) {
                        int32_t swap = faces[3*f+1]; faces[3*f+1] = faces[3*f+2]; faces[3*f+2] = swap;
                    }
                    if (CloseupSurface_audit_proposal_metric(verts,nv,uv,frame,faces,nf,face_patch,1,&passed,&metric)) fails++;
                    int expected = trial == 1 || trial == 3;
                    if (passed != expected || metric.metric_qualified != expected) fails++;
                    fprintf(stderr,"  Provisional metric trial %d: qualified=%d metric=%zu orientation=%zu\n",
                            trial,metric.metric_qualified,metric.metric_failures,metric.orientation_failures);
                }
            }
            free(grid_uv); free(verts); free(frame); free(normals); free(faces); free(nodes); free(uv); free(face_patch);
        }
test_done:
        free(patch); free(original); free(original_state); CloseupGrid_free(&g);
    }
    fprintf(stderr,"  Provisional Cartesian infill: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails != 0;
}
