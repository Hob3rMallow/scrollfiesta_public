/* An endpoint proposal may begin with a few inverted faces. Squared distance
 * to positive, bounded singular values stays finite there, unlike a barrier.
 * The gradient is 2*(J-projection(J)); the final independent certificate,
 * not the optimizer's objective, decides whether this map can be retained. */
static double amb_project(const double j[4], double residual[4])
{
    double ct = j[0]+j[3], sn = j[2]-j[1], length = hypot(ct, sn);
    if (length > 0) { ct /= length; sn /= length; } else { ct = 1; sn = 0; }
    double a = ct*j[0]+sn*j[2], c = -sn*j[1]+ct*j[3];
    double b = .5*(ct*j[1]+sn*j[3]-sn*j[0]+ct*j[2]);
    double gap = hypot(a-c, 2*b), hi = .5*(a+c+gap), lo = .5*(a+c-gap);
    double angle = .5*atan2(2*b, a-c), u = cos(angle), v = sin(angle);
    hi = fmin(1.0998, fmax(.9002, hi)); lo = fmin(1.0998, fmax(.9002, lo));
    double p = hi*u*u+lo*v*v, q = (hi-lo)*u*v, r = hi*v*v+lo*u*u;
    residual[0] = j[0]-(ct*p-sn*q); residual[1] = j[1]-(ct*q-sn*r);
    residual[2] = j[2]-(sn*p+ct*q); residual[3] = j[3]-(sn*q+ct*r);
    double loss = 0; for (int k = 0; k < 4; k++) loss += residual[k]*residual[k];
    return loss;
}

typedef struct AmrBand {
    AsmField *field;
    const size_t *vertices, *faces;
    const int32_t *index;
    const double *seed, *scale;
    size_t n, nf;
    double *contact_gradient, contact_weight;
    AsmContactStats contacts;
    const uint8_t *seams;
    double seam_weight;
} AmrBand;

/* A source path can start with valid seam fits but overlapping interiors.
 * Carry both physical witness directions while the region changes shape. */
static double amb_seam_loss(const AmrBand *s, double *gradient)
{
    if (!s->seams || !(s->seam_weight > 0)) return 0;
    const AsmField *f = s->field; double loss = 0;
    for (size_t i = 0; i < f->n_seams; i++) if (s->seams[i] && !f->seams[i].original_cut) {
        const AsmFieldSeam *seam = f->seams+i; if (!seam->complete) return INFINITY;
        for (size_t k = 0; k < seam->count; k++) for (int d = 0; d < 2; d++) {
            const AsmWitness *w = f->witness+seam->first+k;
            if (!w->valid || !(w->weight > 0)) return INFINITY;
            double r[2]; AsmField_residual(w,d,f->uv,r);
            double length = hypot(r[0],r[1]), excess = fmax(0,length-1.5), weight = .5*s->seam_weight*w->weight;
            if (!(excess > 0)) continue;
            loss += weight*excess*excess;
            for (int q = 0; q < 4; q++) {
                int32_t at = s->index[w->vertex[d][q]]; if (at < 0) continue;
                for (int axis = 0; axis < 2; axis++) gradient[2*(size_t)at+axis] +=
                    2*weight*excess*r[axis]/length*w->coefficient[d][q]/s->scale[at];
            }
        }
    }
    return loss;
}

static int amb_selected_pass(const AsmField *f, const uint8_t *seams)
{
    if (seams) for (size_t i = 0; i < f->n_seams; i++) if (seams[i]) {
        AsmSeamMeasure m; AsmField_seam(f,i,f->uv,&m); if (!m.pass) return 0;
    }
    return 1;
}

static int amb_moving(const uint8_t *core, size_t chart, size_t c)
{ return core ? core[c] != 0 : chart == c; }

static int amb_metrics_pass(const AsmField *f, size_t chart, const uint8_t *core)
{
    for (size_t c = 0; c < f->run->n_charts; c++) if (amb_moving(core,chart,c)) {
        if (f->face_offset[c] == SIZE_MAX) return 0;
        AmrMetric m = amr_metric(f,c,f->uv); if (!amr_qualified(&m)) return 0;
    }
    return 1;
}

static double amb_cut_loss(const AmrBand *s, double *gradient)
{
    const AsmField *f = s->field; double loss = 0;
    for (size_t i = 0; i < f->n_seams; i++) {
        const AsmFieldSeam *seam = f->seams+i;
        if (!seam->original_cut || !seam->complete) continue;
        for (size_t k = 0; k < seam->count; k++) {
            const AsmWitness *w = f->witness+seam->first+k; double r[2]; AsmField_residual(w,0,f->uv,r);
            double length = hypot(r[0],r[1]), excess = fmax(0,length-1.9);
            if (!(excess > 0)) continue;
            loss += 100*excess*excess;
            for (int v = 0; v < 4; v++) {
                int32_t at = s->index[w->vertex[0][v]]; if (at < 0) continue;
                for (int d = 0; d < 2; d++) gradient[2*(size_t)at+d] +=
                    200*excess*r[d]/length*w->coefficient[0][v]/s->scale[at];
            }
        }
    }
    return loss;
}

static int amb_cuts_pass(const AsmField *f)
{
    for (size_t i = 0; i < f->n_seams; i++) if (f->seams[i].original_cut) {
        AsmSeamMeasure m; AsmField_seam(f,i,f->uv,&m); if (!m.pass) return 0;
    }
    return 1;
}

static double amb_loss(AmrBand *s, const double *x, double *gradient)
{
    AsmField *f = s->field; double loss = 0;
    memset(gradient, 0, 2*s->n*sizeof(double));
    for (size_t i = 0; i < s->n; i++) for (int d = 0; d < 2; d++)
        f->uv[2*s->vertices[i]+d] = s->seed[2*i+d]+x[2*i+d]/s->scale[i];
    for (size_t k = 0; k < s->nf; k++) {
        size_t t = s->faces[k]; double j[4], lo, hi, det, residual[4], basis[6];
        AsmField_face(f, t, f->uv, j, &lo, &hi, &det); loss += amb_project(j, residual);
        apr_gradient_basis(f->metric+t, basis);
        for (int v = 0; v < 3; v++) {
            int32_t i = s->index[f->faces[3*t+v]]; if (i < 0) continue;
            gradient[2*(size_t)i] += 2*(residual[0]*basis[2*v]+residual[1]*basis[2*v+1])/s->scale[i];
            gradient[2*(size_t)i+1] += 2*(residual[2]*basis[2*v]+residual[3]*basis[2*v+1])/s->scale[i];
        }
    }
    for (size_t i = 0; i < s->n; i++) for (int d = 0; d < 2; d++) {
        double move = x[2*i+d]/s->scale[i]; loss += 1e-8*move*move;
        gradient[2*i+d] += 2e-8*move/s->scale[i];
    }
    if (s->contact_weight > 0) {
        double overlap_objective=0;
        if (AsmContacts_objective(f->contacts,f->uv,s->contact_gradient,NULL,NULL,&s->contacts,&overlap_objective)) return INFINITY;
        loss += s->contact_weight*overlap_objective;
        for (size_t i = 0; i < s->n; i++) for (int d = 0; d < 2; d++)
            gradient[2*i+d] += s->contact_weight*s->contact_gradient[2*s->vertices[i]+d]/s->scale[i];
        loss += amb_cut_loss(s,gradient);
        loss += amb_seam_loss(s,gradient);
    }
    return isfinite(loss) ? loss : INFINITY;
}

static double amb_dot(const double *a, const double *b, size_t n)
{
    double value = 0; for (size_t i = 0; i < n; i++) value += a[i]*b[i]; return value;
}

typedef struct AmbContactSeeds { const AsmField *field; size_t chart; uint8_t *vertices; const uint8_t *core; } AmbContactSeeds;
static void amb_contact_seed(void *context, size_t a, size_t b, double area)
{
    (void)area; AmbContactSeeds *s = context; size_t face[2] = {a,b};
    for (int side = 0; side < 2; side++) if (amb_moving(s->core,s->chart,(size_t)s->field->face_chart[face[side]]))
        for (int k = 0; k < 3; k++) s->vertices[s->field->faces[3*face[side]+k]] = 1;
}

/* precondition=2 also minimizes complete local contact area. Intermediate
 * proposals may violate constraints; the caller certifies the endpoint. */
static int amb_refine_region_budget(AsmField *f, size_t chart, const uint8_t *core, const uint8_t *seams,
                                   double deadline, int precondition, AmrResult *out, const AsmRepairOpts *opts,
                                   int entire_intrinsic)
{
    Arena_T arena = Arena_new();
    SparseFactor_T factor = NULL;
    int accepted = 0, joint = precondition == 2;
    AsmContactStats initial_contacts = {0};
    if (ves_clock_sec() >= deadline || apr_expired(opts)) { out->reason = "local_metric_deadline"; goto done; }
    if (!apr_workspace_check(opts,apr_workspace_estimate(f->nv,f->nf,f->n_witness,f->run->n_charts))) {
        out->reason = "workspace_memory_budget"; goto done;
    }
    if (!joint && amb_metrics_pass(f,chart,core)) { Arena_dispose(&arena); return 1; }
    uint8_t *free_vertex = ARENA_CALLOC(arena, f->nv, 1), *next = ARENA_ALLOC(arena, f->nv);
    if (joint) {
        if (!f->contacts && AsmField_build_contacts(f)) { out->reason = "joint_contact_tree_failed"; goto done; }
        uint8_t *active = ARENA_CALLOC(arena,f->nf,1);
        for (size_t t = 0; t < f->nf; t++) active[t] = (uint8_t)amb_moving(core,chart,(size_t)f->face_chart[t]);
        AsmContacts_active_faces(f->contacts,active);
        AmbContactSeeds seeds = {f,chart,free_vertex,core};
        if (AsmContacts_measure(f->contacts,f->uv,NULL,amb_contact_seed,&seeds,&initial_contacts)) { out->reason = "joint_contact_walk_failed"; goto done; }
        if (!apr_workspace_check(opts,0)) { out->reason = "workspace_memory_budget"; goto done; }
        if (amb_metrics_pass(f,chart,core) && !initial_contacts.pairs && amb_cuts_pass(f) && amb_selected_pass(f,seams)) { accepted = 1; goto done; }
        for (size_t i = 0; i < f->n_seams; i++) if ((f->seams[i].original_cut || (seams && seams[i])) && f->seams[i].complete) {
            AsmSeamMeasure m; AsmField_seam(f,i,f->uv,&m); if (m.pass) continue;
            for (size_t k = 0; k < f->seams[i].count; k++) {
                const AsmWitness *w = f->witness+f->seams[i].first+k;
                for (int d = 0; d < 2; d++) for (int q = 0; q < 4; q++)
                    if (amb_moving(core,chart,(size_t)f->vertex_chart[w->vertex[d][q]])) free_vertex[w->vertex[d][q]] = 1;
            }
        }
    }
    for (size_t c = 0; c < f->run->n_charts; c++) if (amb_moving(core,chart,c)) {
        /* An unregistered intrinsic seed has no fixed surrounding UV frame.
         * Free its complete chart under the byte/time budget; freezing a
         * topological ring would introduce an artificial fitting constraint. */
        if (entire_intrinsic) {
            memset(free_vertex+f->vertex_offset[c],1,f->run->charts[c].nv); continue;
        }
        AmrMetric initial = amr_metric(f,c,f->uv); int normal_band = initial.within10 < .95*initial.area;
        for (size_t k = 0; k < f->run->charts[c].nf; k++) {
            size_t t = f->face_offset[c]+k; double j[4], lo, hi, det; AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
            if (!(det > 0) || !isfinite(lo) || !isfinite(hi) || lo < (normal_band ? .9 : .75) || hi > (normal_band ? 1.1 : 1.25))
                for (int v = 0; v < 3; v++) free_vertex[f->faces[3*t+v]] = 1;
        }
    }
    for (int ring = 0; !entire_intrinsic && ring < 8; ring++) {
        memcpy(next, free_vertex, f->nv);
        for (size_t t = 0; t < f->nf; t++) if (amb_moving(core,chart,(size_t)f->face_chart[t])) {
            const int32_t *tri = f->faces+3*t;
            if (free_vertex[tri[0]] || free_vertex[tri[1]] || free_vertex[tri[2]])
                for (int v = 0; v < 3; v++) next[tri[v]] = 1;
        }
        uint8_t *swap = next; next = free_vertex; free_vertex = swap;
    }
    size_t n = 0, nf = 0;
    for (size_t v = 0; v < f->nv; v++) n += free_vertex[v] != 0;
    for (size_t t = 0; t < f->nf; t++) {
        const int32_t *tri = f->faces+3*t;
        nf += free_vertex[tri[0]] || free_vertex[tri[1]] || free_vertex[tri[2]];
    }
    out->free_vertices = n; out->active_faces = nf;
    if (!n || n > INT_MAX || nf > (size_t)INT_MAX/6) { out->reason = "local_metric_index_range"; goto done; }
    /* L-BFGS history and sparse assembly are actual workspace. Density is
     * not a geometric eligibility criterion. Native factor fill is sampled
     * separately, as in the enclosing admission transaction. */
    size_t estimated = apr_workspace_estimate(f->nv,f->nf,f->n_witness,f->run->n_charts);
    apr_bytes_add(&estimated,n,80*sizeof(double)+sizeof(size_t));
    apr_bytes_add(&estimated,nf,sizeof(size_t));
    if (!apr_workspace_check(opts,estimated)) { out->reason = "workspace_memory_budget"; goto done; }
    size_t *vertices = ARENA_ALLOC(arena, n*sizeof(size_t)), *faces = ARENA_ALLOC(arena, nf*sizeof(size_t));
    int32_t *index = ARENA_ALLOC(arena, f->nv*sizeof(int32_t));
    for (size_t v = 0; v < f->nv; v++) index[v] = -1;
    size_t at = 0;
    for (size_t v = 0; v < f->nv; v++) if (free_vertex[v]) { index[v] = (int32_t)at; vertices[at++] = v; }
    double *scale = ARENA_CALLOC(arena, n, sizeof(double)), *seed = ARENA_ALLOC(arena, 2*n*sizeof(double));
    at = 0;
    for (size_t t = 0; t < f->nf; t++) {
        const int32_t *tri = f->faces+3*t;
        if (!free_vertex[tri[0]] && !free_vertex[tri[1]] && !free_vertex[tri[2]]) continue;
        faces[at++] = t; double basis[6]; apr_gradient_basis(f->metric+t, basis);
        for (int v = 0; v < 3; v++) if (index[tri[v]] >= 0)
            scale[index[tri[v]]] += basis[2*v]*basis[2*v]+basis[2*v+1]*basis[2*v+1];
    }
    for (size_t i = 0; i < n; i++) { scale[i] = sqrt(fmax(scale[i], 1e-12)); memcpy(seed+2*i, f->uv+2*vertices[i], 2*sizeof(double)); }
    size_t count = 2*n; enum {HISTORY=15};
    double *x = ARENA_CALLOC(arena, count, sizeof(double)), *g = ARENA_ALLOC(arena, count*sizeof(double));
    double *tx = ARENA_ALLOC(arena, count*sizeof(double)), *tg = ARENA_ALLOC(arena, count*sizeof(double));
    double *direction = ARENA_ALLOC(arena, count*sizeof(double));
    double *ss = ARENA_ALLOC(arena, HISTORY*count*sizeof(double)), *yy = ARENA_ALLOC(arena, HISTORY*count*sizeof(double));
    double rho[HISTORY], coefficients[HISTORY]; int history = 0;
    AmrBand ctx = {f,vertices,faces,index,seed,scale,n,nf};
    if (joint) { ctx.contact_weight = 10; ctx.contact_gradient = ARENA_ALLOC(arena,2*f->nv*sizeof(double)); }
    ctx.seams = seams; ctx.seam_weight = seams ? 100 : 0;
    double value = amb_loss(&ctx, x, g);
    double *rhs = NULL, *answer = NULL;
    if (precondition && ves_clock_sec() < deadline) {
        double started = ves_clock_sec(); size_t capacity = 6*nf+n+(seams ? 20*f->n_witness : 0); AprCoo matrix = {0};
        matrix.row = ARENA_ALLOC(arena,capacity*sizeof(int)); matrix.col = ARENA_ALLOC(arena,capacity*sizeof(int));
        matrix.value = ARENA_ALLOC(arena,capacity*sizeof(double));
        for (size_t k = 0; k < nf; k++) {
            size_t t = faces[k]; double basis[6]; apr_gradient_basis(f->metric+t,basis);
            for (int a = 0; a < 3; a++) for (int b = 0; b <= a; b++) {
                int i = index[f->faces[3*t+a]], j = index[f->faces[3*t+b]];
                if (i >= 0 && j >= 0) apr_coo(&matrix,i,j,
                    2*(basis[2*a]*basis[2*b]+basis[2*a+1]*basis[2*b+1])/(scale[i]*scale[j]));
            }
        }
        if (seams) for (size_t si = 0; si < f->n_seams; si++) if (seams[si] && f->seams[si].complete) {
            const AsmFieldSeam *seam = f->seams+si;
            for (size_t k = 0; k < seam->count; k++) for (int d = 0; d < 2; d++) {
                const AsmWitness *w = f->witness+seam->first+k;
                for (int a = 0; a < 4; a++) for (int b = 0; b <= a; b++) {
                    int i = index[w->vertex[d][a]], j = index[w->vertex[d][b]];
                    if (i >= 0 && j >= 0) apr_coo(&matrix,i,j,ctx.seam_weight*w->weight*w->coefficient[d][a]*w->coefficient[d][b]/(scale[i]*scale[j]));
                }
            }
        }
        for (size_t i = 0; i < n; i++) apr_coo(&matrix,(int)i,(int)i,2e-8/(scale[i]*scale[i]));
        int rc = Sparse_factor_spd((int)n,matrix.n,matrix.row,matrix.col,matrix.value,&factor);
        out->factor_seconds = ves_clock_sec()-started;
        if (rc) { out->reason = "local_metric_factor_failed"; goto done; }
        if (!apr_workspace_check(opts,0)) { out->reason = "workspace_memory_budget"; goto done; }
        rhs = ARENA_ALLOC(arena,count*sizeof(double)); answer = ARENA_ALLOC(arena,count*sizeof(double));
    }
    out->reason = "local_metric_stalled";
    while (isfinite(value) && ves_clock_sec() < deadline) {
        if (!apr_workspace_check(opts,0)) { out->reason = "workspace_memory_budget"; break; }
        memcpy(direction, g, count*sizeof(double));
        for (int h = history-1; h >= 0; h--) {
            coefficients[h] = rho[h]*amb_dot(ss+(size_t)h*count, direction, count);
            for (size_t i = 0; i < count; i++) direction[i] -= coefficients[h]*yy[(size_t)h*count+i];
        }
        if (factor) {
            double started = ves_clock_sec();
            for (size_t i = 0; i < n; i++) for (int d = 0; d < 2; d++) rhs[(size_t)d*n+i] = direction[2*i+d];
            int rc = Sparse_factor_solve_multi(factor,rhs,answer,2); out->solve_seconds += ves_clock_sec()-started;
            if (rc) { out->reason = "local_metric_solve_failed"; break; }
            for (size_t i = 0; i < n; i++) for (int d = 0; d < 2; d++) direction[2*i+d] = answer[(size_t)d*n+i];
        } else if (history) {
            double norm = amb_dot(yy+(size_t)(history-1)*count, yy+(size_t)(history-1)*count, count);
            double multiplier = 1/(rho[history-1]*norm);
            for (size_t i = 0; i < count; i++) direction[i] *= multiplier;
        }
        for (int h = 0; h < history; h++) {
            double beta = rho[h]*amb_dot(yy+(size_t)h*count, direction, count);
            for (size_t i = 0; i < count; i++) direction[i] += ss[(size_t)h*count+i]*(coefficients[h]-beta);
        }
        for (size_t i = 0; i < count; i++) direction[i] = -direction[i];
        double slope = amb_dot(direction, g, count);
        if (!(slope < 0) || !isfinite(slope)) {
            history = 0; for (size_t i = 0; i < count; i++) direction[i] = -g[i]; slope = -amb_dot(g, g, count);
        }
        if (!(slope < 0) || !isfinite(slope)) break;
        double after = INFINITY; int found = 0;
        for (int line = 0; line < 40 && ves_clock_sec() < deadline; line++) {
            double alpha = ldexp(1., -line);
            for (size_t i = 0; i < count; i++) tx[i] = x[i]+alpha*direction[i];
            after = amb_loss(&ctx, tx, tg);
            if (after < value && after <= value+1e-4*alpha*slope) { found = 1; break; }
        }
        if (!found) break;
        out->refine_steps++;
        if (amb_metrics_pass(f,chart,core) && (!joint || (!ctx.contacts.pairs && amb_cuts_pass(f) && amb_selected_pass(f,seams)))) { accepted = 1; break; }
        double sy = 0, snorm = 0, ynorm = 0;
        for (size_t i = 0; i < count; i++) {
            double a = tx[i]-x[i], b = tg[i]-g[i]; sy += a*b; snorm += a*a; ynorm += b*b;
        }
        if (isfinite(sy) && sy > 1e-12*sqrt(snorm*ynorm)) {
            if (history == HISTORY) {
                memmove(ss, ss+count, (HISTORY-1)*count*sizeof(double)); memmove(yy, yy+count, (HISTORY-1)*count*sizeof(double));
                memmove(rho, rho+1, (HISTORY-1)*sizeof(double)); history--;
            }
            for (size_t i = 0; i < count; i++) { ss[(size_t)history*count+i] = tx[i]-x[i]; yy[(size_t)history*count+i] = tg[i]-g[i]; }
            rho[history++] = 1/sy;
        }
        memcpy(x, tx, count*sizeof(double)); memcpy(g, tg, count*sizeof(double)); value = after;
    }
    if (!accepted) amb_loss(&ctx, x, g);
    if (ves_clock_sec() >= deadline) { accepted = 0; out->reason = "local_metric_deadline"; }
done:
    if (accepted) out->reason = "local_metric_qualified";
    if (joint && f->contacts) AsmContacts_active_faces(f->contacts,NULL);
    Sparse_factor_free(&factor); Arena_dispose(&arena); return accepted;
}

static int amb_refine_region(AsmField *f, size_t chart, const uint8_t *core, const uint8_t *seams,
                             double deadline, int precondition, AmrResult *out)
{
    AsmRepairOpts opts; AsmRepair_defaults(&opts); opts.compact = 1; opts.deadline = deadline;
    AprWorkspace work = {0}; work.measured = apr_process_bytes(&work.entry); opts.workspace = &work;
    return amb_refine_region_budget(f,chart,core,seams,deadline,precondition,out,&opts,0);
}

static int amb_refine(AsmField *f, size_t chart, double deadline, int precondition, AmrResult *out)
{ return amb_refine_region(f,chart,NULL,NULL,deadline,precondition,out); }
