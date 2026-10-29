/* A cone is impossible while closed, but has an exact development after one
 * radial cut. Verify the construction against that analytic development,
 * source identity, frozen witness branches, contacts, and serialized replay. */

static int acut_clearance_selftest(void)
{
    int failures = 0;
    /* Vertex polish corrects a sub-voxel excess over the contact band; its
     * physical steps share the repair deadline. Cases 4 and 5 (mirrored). */
    for (int test = 4; test < 6; test++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[] = {0,0,0, 1,0,0, 0,1,0, 5,0,0, 6,0,0, 5,1,0};
        float uv[] = {0,0, 1,0, 0,1, .85f,0, 1.85f,0, .85f,1};
        /* After scaling, diameter (2-sqrt(2))*.108*64 = 4.049 vox. */
        uv[6]=uv[10]=.892f;uv[8]=1.892f;
        int32_t faces[] = {0,1,2,3,4,5}; double placed[12];
        AsmChart ch = {0}; ch.id = 0; ch.parent = -1; ch.nv = 6; ch.nf = 2;
        ch.xyz = xyz; ch.uv = uv; ch.faces = faces; ch.placed_uv = placed;
        ch.placed = 1; ch.placement_state = ASM_PLACE_ROOT;
        if (test%2) ch.flags |= ASM_CHART_MIRROR;
        for (int v = 0; v < 6; v++) {
            placed[2*v] = (test%2?1e9:0)+(test%2?-1:1)*.5*xyz[3*v];
            placed[2*v+1] = (test%2?-1e9:0)+.5*xyz[3*v+1];
        }
        run.charts = &ch; run.n_charts = 1;
        apr_test_scale(&run,64); /* The obstruction must exceed the 4-voxel band. */
        AsmField *f = AsmField_new(arena,&run); AmrResult r = {0}; AsmRepairOpts opts; AsmRepair_defaults(&opts);
        opts.patch_faces = opts.patch_field_faces = 2; opts.deadline = ves_clock_sec()+5;
        int okay = 0;
        if (f) {
            for (int v = 0; v < 6; v++) {
                f->uv[2*v] = (test%2?1e9:0)+(test%2?-1:1)*uv[2*v];
                f->uv[2*v+1] = (test%2?-1e9:0)+uv[2*v+1];
            }
            r.clearance_steps = amr_vertex_polish(f,0,opts.deadline); r.reason = "vertex_polish";
            AsmContactStats contacts = {0}; r.after = amr_metric(f,0,f->uv);
            okay = r.clearance_steps && amr_qualified(&r.after) && !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&contacts) && !contacts.pairs;
        }
        fprintf(stderr,"  native metric endpoint local clearance %d: %s, %zu steps, %s\n",test,r.reason?r.reason:"field_failed",r.clearance_steps,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    return failures;
}

static int acut_boundary_selftest(void)
{
    enum {N=9, NV=N*N, NF=2*(N-1)*(N-1)}; int failures = 0;
    for (int outer = 0; outer < 3; outer++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}, trial = {0}; run.arena = arena;
        float xyz[3*NV], uv[2*NV]; int32_t faces[3*NF], vid[NV]; size_t nf = 0;
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
            int v = y*N+x; vid[v] = v; xyz[3*v] = uv[2*v] = (float)x; xyz[3*v+1] = uv[2*v+1] = (float)y;
            xyz[3*v+2] = x == 4 && y == 4 ? 4.f : 0.f;
            if (x+1 < N && y+1 < N && !(x == 5 && y == 4)) {
                int32_t tri[] = {v,v+1,v+N+1,v,v+N+1,v+N}; memcpy(faces+3*nf,tri,sizeof tri); nf += 2;
            }
        }
        AsmChart ch = {0}; ch.id = 0; ch.parent = -1; ch.nv = NV; ch.nf = nf;
        ch.xyz = xyz; ch.uv = uv; ch.vid = vid; ch.faces = faces; ch.placed = 1; ch.placement_state = ASM_PLACE_ROOT;
        AsmChart charts[2] = {ch,ch}; charts[1].id = 1; run.charts = charts; run.n_charts = 2;
        AsmRelation relation = {0}; AsmCorr corr[8] = {{0}};
        const int32_t sampled[8] = {0,4,8,4*N+8,8*N+8,8*N+4,8*N,4*N};
        relation.a = 0; relation.b = 1; relation.flags = ASM_REL_WEAK;
        relation.continuity = ASM_CONT_SOURCE | ASM_CONT_REQUIRED; relation.n_corr = relation.corr_count = 8;
        run.rels = &relation; run.n_rels = 1; run.corr = corr; run.n_corr = 8;
        for (int k = 0; k < 8; k++) {
            corr[k].va = corr[k].vb = sampled[k]; corr[k].valid = 3;
            corr[k].trim_a.valid = corr[k].trim_b.valid = 1;
            for (size_t face = 0; face < nf; face++) for (int corner = 0; corner < 3; corner++)
                if (faces[3*face+corner] == sampled[k]) corr[k].trim_a.face = corr[k].trim_b.face = (int32_t)face;
        }
        AcutResult r;
        int okay = acut_make_mode(&run,0,4*N+4,arena,&trial,&r,outer) && r.source_preserved &&
            r.boundary_components == 2 && fabs(r.target_perimeter-(outer == 2 ? 36 : outer?32:4)) < 1e-12;
        AsmField *f = okay ? AsmField_new_local(arena,&trial,NULL,NULL) : NULL;
        okay = okay && f && f->n_seams == r.edges+1;
        if (f) for (size_t i = 0; i < f->n_seams; i++) { AsmSeamMeasure seam; AsmField_seam(f,i,f->uv,&seam); okay = okay && seam.pass; }
        if (f) {
            CSR_T graph = CSR_from_faces(arena,trial.charts->faces,nf,trial.charts->nv);
            int64_t euler = (int64_t)trial.charts->nv-CSR_nnz(graph)/2+(int64_t)nf;
            okay = okay && euler == (outer == 2 ? 1 : 0);
        }
        fprintf(stderr,"  native source cut boundary selection %d: %s, %zu edges, %s\n",outer,r.reason,r.edges,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    return failures;
}

static int acut_density_selftest(void)
{
    enum {N=8}; int failures = 0; double coarse_length = 0;
    for (int test = 0; test < 4; test++) {
        int rings = test ? 65 : 1; size_t nv = 1+(size_t)N*rings, nf = (size_t)N*(2*rings-1);
        Arena_T arena = Arena_new(); AsmRun run = {0}, trial = {0}; run.arena = arena;
        float *xyz = ARENA_CALLOC(arena,3*nv,sizeof(float)), *uv = ARENA_CALLOC(arena,2*nv,sizeof(float));
        int32_t *vid = ARENA_ALLOC(arena,nv*sizeof(int32_t)), *faces = ARENA_ALLOC(arena,3*nf*sizeof(int32_t));
        xyz[2] = 1; vid[0] = 0; size_t face = 0;
        for (int ring = 1; ring <= rings; ring++) for (int k = 0; k < N; k++) {
            int32_t v = 1+(ring-1)*N+k; double radius = (double)ring/rings, angle = 2*acos(-1.)*k/N;
            vid[v] = v; xyz[3*v] = uv[2*v] = (float)(.25*radius*cos(angle));
            xyz[3*v+1] = uv[2*v+1] = (float)(.25*radius*sin(angle)); xyz[3*v+2] = (float)(1-radius);
            if (ring == 1) {
                faces[3*face] = 0; faces[3*face+1] = v; faces[3*face++ +2] = 1+(k+1)%N;
            } else {
                int32_t a = v-N, b = 1+(ring-2)*N+(k+1)%N, c = 1+(ring-1)*N+(k+1)%N;
                int32_t tri[] = {a,v,c,a,c,b}; memcpy(faces+3*face,tri,sizeof tri); face += 2;
            }
        }
        AsmChart ch = {0}; ch.id = 0; ch.parent = -1; ch.nv = nv; ch.nf = nf;
        ch.xyz = xyz; ch.uv = uv; ch.vid = vid; ch.faces = faces; ch.placed = 1; ch.placement_state = ASM_PLACE_ROOT;
        run.charts = &ch; run.n_charts = 1;
        AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1; opts.deadline = ves_clock_sec()+5;
        AprWorkspace work = {0}; opts.workspace = &work;
        if (test == 2) opts.patch_bytes = 1;
        if (test == 3) opts.deadline = ves_clock_sec()-1;
        AcutResult r; int made = acut_make_mode_budget(&run,0,SIZE_MAX,arena,&trial,&r,0,&opts);
        int okay = test < 2 ? made : !made && !strcmp(r.reason,test == 2 ? "workspace_memory_budget" : "cut_deadline");
        if (made) {
            okay = okay && face == nf && r.edges == (size_t)rings && r.added_vertices == (size_t)rings &&
                r.source_preserved && r.orientation_before == 1 && !r.orientation_after;
            if (!test) coarse_length = r.length;
            else okay = okay && fabs(r.length-coarse_length) < 1e-6;
            AsmChart *next = trial.charts; int first = -1;
            for (int k = 0; k < N; k++) {
                int incoming = 0;
                for (int j = 0; j < N; j++) incoming += next->faces[3*j+2] == next->faces[3*k+1];
                if (!incoming) first = k;
            }
            AsmField *f = AsmField_new_local(arena,&trial,NULL,NULL); okay = okay && f && first >= 0;
            if (f && first >= 0) {
                double length = sqrt(1.0625), edge = .5*sin(acos(-1.)/N);
                double wedge = 2*asin(edge/(2*length));
                for (size_t i = 0; i < nf; i++) for (int corner = 0; corner < 3; corner++) {
                    int32_t v = next->faces[3*i+corner], original = next->vid[v];
                    if (!original) { f->uv[2*v] = f->uv[2*v+1] = 0; continue; }
                    int angular = (original-1)%N, ring = (original-1)/N+1;
                    int sector = i < N ? (int)i : (int)((i-N)/2%N);
                    int turn = (angular-first+N)%N;
                    if (!turn && (sector+1)%N == first) turn = N;
                    double radius = length*ring/rings;
                    f->uv[2*v] = radius*cos(wedge*turn); f->uv[2*v+1] = radius*sin(wedge*turn);
                }
                AmrMetric metric = amr_metric(f,0,f->uv); okay = okay && amr_qualified(&metric);
                for (size_t i = 0; i < f->n_seams; i++) {
                    AsmSeamMeasure seam; AsmField_seam(f,i,f->uv,&seam); okay = okay && seam.pass;
                }
                AsmContactStats contacts = {0};
                okay = okay && !AsmField_build_contacts(f) && !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&contacts) && !contacts.pairs;
                CSR_T graph = CSR_from_faces(arena,next->faces,nf,next->nv);
                okay = okay && (int64_t)next->nv-CSR_nnz(graph)/2+(int64_t)nf == 1;
            }
        }
        fprintf(stderr,"  native source cut density %d: %s, %zu edges, %.9g source length, %s\n",test,r.reason,r.edges,r.length,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    return failures;
}

static int acut_forest_selftest(void)
{
    enum {N=13,NV=N*N,NF=2*(N-1)*(N-1)};
    Arena_T arena = Arena_new(); AsmRun run = {0}, trial = {0}; run.arena = arena;
    float xyz[3*NV], uv[2*NV]; int32_t faces[3*NF], vid[NV]; size_t nf = 0;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        int v = y*N+x; vid[v] = v; xyz[3*v] = uv[2*v] = (float)x; xyz[3*v+1] = uv[2*v+1] = (float)y;
        xyz[3*v+2] = (x == 4 && y == 4) || (x == 8 && y == 8) ? 8.f : 0.f;
        if (x+1 < N && y+1 < N) { int32_t tri[] = {v,v+1,v+N+1,v,v+N+1,v+N}; memcpy(faces+3*nf,tri,sizeof tri); nf += 2; }
    }
    AsmChart ch = {0}; ch.id = 0; ch.parent = -1; ch.nv = NV; ch.nf = nf;
    ch.xyz = xyz; ch.uv = uv; ch.vid = vid; ch.faces = faces; ch.placed = 1; ch.placement_state = ASM_PLACE_ROOT;
    run.charts = &ch; run.n_charts = 1; AcutResult r;
    int okay = acut_make(&run,0,SIZE_MAX,arena,&trial,&r) && r.source_preserved && r.orientation_before >= 2 && !r.orientation_after;
    if (okay) {
        AsmChart *next = trial.charts; CSR_T graph = CSR_from_faces(arena,next->faces,nf,next->nv);
        size_t *root = ARENA_ALLOC(arena,nf*sizeof(size_t)); for (size_t i = 0; i < nf; i++) root[i] = i;
        size_t *edge_face = ARENA_ALLOC(arena,(size_t)CSR_nnz(graph)*sizeof(size_t));
        for (size_t k = 0; k < (size_t)CSR_nnz(graph); k++) edge_face[k] = SIZE_MAX;
        const int32_t *offset = CSR_offset(graph), *target = CSR_target(graph);
        for (size_t k = 0; k < 3*nf; k++) edge_face[acut_slot(offset,target,next->faces[k],next->faces[acut_next(k)])] = k/3;
        for (size_t k = 0; k < 3*nf; k++) {
            int32_t slot = acut_slot(offset,target,next->faces[acut_next(k)],next->faces[k]);
            if (edge_face[slot] != SIZE_MAX) root[apr_root(root,k/3)] = apr_root(root,edge_face[slot]);
        }
        for (size_t i = 0; i < nf; i++) okay = okay && apr_root(root,i) == apr_root(root,0);
        okay = okay && (int64_t)next->nv-CSR_nnz(graph)/2+(int64_t)nf == 1;
    }
    fprintf(stderr,"  native multiple-obstruction connected cut forest: %s, %zu -> %zu obstructions, %zu edges, %s\n",
        r.reason,r.orientation_before,r.orientation_after,r.edges,okay?"ok":"FAIL");
    Arena_dispose(&arena); return !okay;
}

static int acut_selftest(void)
{
    enum {N=8}; int failures = acut_clearance_selftest()+acut_boundary_selftest()+acut_density_selftest()+acut_forest_selftest();
    for (int test = 0; test < 9; test++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}, trial = {0}; run.arena = arena; run.continuity_ready = 1;
        float xyz[3*(N+1)] = {0}, uv[2*(N+1)] = {0}; int32_t vid[N+1], faces[3*(N+1)] = {0};
        xyz[2] = test == 2 ? 0.f : 1.f;
        for (int v = 0; v <= N; v++) vid[v] = v;
        for (int k = 0; k < N; k++) {
            double angle = 2*acos(-1.)*k/N;
            xyz[3*(k+1)] = uv[2*(k+1)] = (float)(.25*cos(angle));
            xyz[3*(k+1)+1] = uv[2*(k+1)+1] = (float)(.25*sin(angle));
            faces[3*k] = 0; faces[3*k+1] = k+1; faces[3*k+2] = (k+1)%N+1;
        }
        AsmChart charts[2] = {{0}}; AsmChart *ch = charts;
        ch->id = 0; ch->parent = -1; ch->nv = N+1; ch->nf = N;
        ch->xyz = xyz; ch->uv = uv; ch->vid = vid; ch->faces = faces; ch->placed = 1; ch->placement_state = ASM_PLACE_ROOT;
        if (test == 1) { ch->flags |= ASM_CHART_MIRROR; ch->pose_x = 1e9; ch->pose_y = -1e9; }
        if (test == 5) { memcpy(faces+3*N,faces,3*sizeof(int32_t)); ch->nf++; }
        if (test == 6) ch->nf = (size_t)INT32_MAX/6+1;
        charts[1].id = 1; run.charts = charts; run.n_charts = 2;
        AsmRelation relation = {0}; AsmCorr corr[2*N] = {{0}};
        if (test == 4) {
            relation.a = 0; relation.b = 1; relation.flags = ASM_REL_WEAK; relation.continuity = ASM_CONT_SOURCE;
            relation.n_corr = relation.corr_count = 2*N; run.rels = &relation; run.n_rels = 1; run.corr = corr; run.n_corr = 2*N;
            for (int k = 0; k < 2*N; k++) {
                corr[k].va = faces[3*(k/2)+1+k%2]; corr[k].trim_a.valid = 1; corr[k].trim_a.face = k/2;
            }
        }
        if (test == 8) {
            charts[1] = *ch; charts[1].id = 1;
            relation.a = 0; relation.b = 1; relation.flags = ASM_REL_WEAK; relation.continuity = ASM_CONT_SOURCE | ASM_CONT_REQUIRED;
            relation.n_corr = relation.corr_count = N; run.rels = &relation; run.n_rels = 1; run.corr = corr; run.n_corr = N;
            for (int k = 0; k < N; k++) {
                corr[k].va = corr[k].vb = k+1; corr[k].valid = 3;
                corr[k].trim_a.valid = corr[k].trim_b.valid = 1; corr[k].trim_a.face = corr[k].trim_b.face = k;
            }
        }
        AsmChart saved_charts[2]; float saved_xyz[sizeof xyz/sizeof *xyz]; int32_t saved_faces[sizeof faces/sizeof *faces];
        memcpy(saved_charts,charts,sizeof charts); memcpy(saved_xyz,xyz,sizeof xyz); memcpy(saved_faces,faces,sizeof faces);
        AcutResult result; int made = acut_make(&run,0,test == 3 ? N+1 : 0,arena,&trial,&result);
        int expected = test == 0 || test == 1 || test == 4 || test == 7 || test == 8;
        int okay = made == expected && !memcmp(saved_charts,charts,sizeof charts) &&
            !memcmp(saved_xyz,xyz,sizeof xyz) && !memcmp(saved_faces,faces,sizeof faces);
        if (made) {
            AsmChart *next = trial.charts; okay = okay && result.source_preserved && result.edges == 1 &&
                result.added_vertices == 1 && result.before_obstructions == 1 && !result.after_obstructions;
            if (test == 4) {
                int changed = 0;
                for (int k = 0; k < 2*N; k++) {
                    okay = okay && trial.corr[k].va == next->faces[3*(k/2)+1+k%2];
                    changed += trial.corr[k].va != corr[k].va;
                }
                okay = okay && changed > 0 && !memcmp(&relation,trial.rels,sizeof relation);
            }
            AsmField *f = AsmField_new_local(arena,&trial,NULL,NULL);
            okay = okay && f;
            if (f) {
                if (test == 8) {
                    AsmSeamMeasure original; AsmField_seam(f,0,f->uv,&original);
                    okay = okay && f->seams[0].complete && !original.invalid && original.observations == N &&
                        original.mass > 0 && original.support == original.mass && original.rms == 0;
                }
                AsmSeamMeasure seam; AsmField_seam(f,f->n_seams-1,f->uv,&seam); okay = okay && seam.pass;
                /* Follow the opened outer chain; preserve each source angle. */
                int first = -1;
                for (int k = 0; k < N; k++) {
                    int incoming = 0;
                    for (int j = 0; j < N; j++) incoming += next->faces[3*j+2] == next->faces[3*k+1];
                    if (!incoming) first = k;
                }
                okay = okay && first >= 0;
                if (first >= 0) {
                    double angle = 0, sign = test == 1 ? -1 : 1;
                    f->uv[0] = ch->pose_x; f->uv[1] = ch->pose_y;
                    for (int k = 0; k < N; k++) {
                        int face = (first+k)%N; size_t a = (size_t)next->faces[3*face+1], b = (size_t)next->faces[3*face+2];
                        double la = acut_length(next,0,a), lb = acut_length(next,0,b), edge = acut_length(next,a,b);
                        f->uv[2*a] = ch->pose_x+la*cos(angle); f->uv[2*a+1] = ch->pose_y+sign*la*sin(angle);
                        angle += acos(fmax(-1,fmin(1,(la*la+lb*lb-edge*edge)/(2*la*lb))));
                        f->uv[2*b] = ch->pose_x+lb*cos(angle); f->uv[2*b+1] = ch->pose_y+sign*lb*sin(angle);
                    }
                    AmrMetric metric = amr_metric(f,0,f->uv); AsmField_seam(f,f->n_seams-1,f->uv,&seam);
                    AsmContactStats contact = {0};
                    okay = okay && amr_qualified(&metric) && seam.pass && !AsmField_build_contacts(f) &&
                        !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&contact) && !contact.pairs;
                    if (test == 7) {
                        AsmRun replay = {0}; replay.arena = arena;
                        const char *path = "asm_metric_cut_selftest.asr";
                        okay = okay && !AsmField_commit(f,&trial) && !AsmStore_write_run(path,&trial) && !AsmStore_read_run(path,&replay);
                        AsmField *loaded = replay.n_charts ? AsmField_new_local(arena,&replay,NULL,NULL) : NULL;
                        okay = okay && loaded && acut_source_same(ch,replay.charts);
                        if (loaded) { AsmField_seam(loaded,loaded->n_seams-1,loaded->uv,&seam); metric = amr_metric(loaded,0,loaded->uv); okay = okay && seam.pass && amr_qualified(&metric); }
                        if (loaded) {
                            int same = acut_old_seams_same(&trial,&replay,0,0);
                            const AsmCorr *c = replay.corr+replay.rels[0].corr_first;
                            if (c->va == c->vb) c++;
                            replay.charts[0].placed_uv[2*(size_t)c->vb] += .125;
                            int changed = !acut_old_seams_same(&trial,&replay,0,0);
                            fprintf(stderr,"  native pre-cut witness preservation and perturbation refusal: %s\n",same && changed?"ok":"FAIL");
                            okay = okay && same && changed;
                        }
                        remove(path);
                    }
                    if (test == 0) {
                        /* Independent finite differences of the exact cut's
                         * endpoint-gap penalty, including a shared slit tip. */
                        int32_t *index = ARENA_ALLOC(arena,f->nv*sizeof(int32_t));
                        double *scale = ARENA_ALLOC(arena,f->nv*sizeof(double));
                        double *gradient = ARENA_CALLOC(arena,2*f->nv,sizeof(double));
                        for (size_t v = 0; v < f->nv; v++) { index[v] = (int32_t)v; scale[v] = 1; }
                        AmrBand band = {0}; band.field = f; band.index = index; band.scale = scale;
                        f->uv[2*(f->nv-1)] += 3;
                        double value = amb_cut_loss(&band,gradient); int derivative = value > 0;
                        for (size_t k = 0; k < 2*f->nv; k++) {
                            double saved = f->uv[k]; double *scratch = ARENA_CALLOC(arena,2*f->nv,sizeof(double));
                            f->uv[k] = saved+1e-6; double a = amb_cut_loss(&band,scratch);
                            f->uv[k] = saved-1e-6; double b = amb_cut_loss(&band,scratch); f->uv[k] = saved;
                            derivative = derivative && fabs((a-b)/2e-6-gradient[k]) < 1e-5*fmax(1,fabs(gradient[k]));
                        }
                        fprintf(stderr,"  native original-cut penalty gradient: %s\n",derivative?"ok":"FAIL");
                        okay = okay && derivative;
                    }
                }
            }
        }
        fprintf(stderr,"  native source-preserving metric cut %d: %s, %s\n",test,result.reason,okay?"ok":"FAIL");
        failures += !okay; Arena_dispose(&arena);
    }
    return failures;
}
