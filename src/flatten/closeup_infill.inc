/* Cartesian thin-plate continuation. Unlike the cylindrical ribbon solve,
 * this operates on the actual sheet frame, in all three world coordinates.
 * The first pass is anchored only to observations. A second pass closes
 * enclosed gaps in already generated surface without moving its collar or
 * relabelling generated geometry as observed papyrus. */
#define CS_INFILL_MAX_NODES 65536u

typedef struct CsInfillSystem {
    const CloseupGrid *g;
    const double *uv;
    const int32_t *index;
    const int32_t *owner;          /* with owner: index is fixed for every hole; a node is this hole's when owner says */
    int32_t region;
    int32_t component;
    int fill_collar;
    int n, nt, cap, dimensions;
    int *rows, *cols;
    double *values, *rhs;
    double origin[5];
} CsInfillSystem;

/* The node's unknown in this hole's system, or -1. */
static int32_t cs_infill_idx(const CsInfillSystem *s, size_t at)
{
    if (s->owner != NULL && s->owner[at] != s->region) return -1;
    return s->index[at];
}

/* A hole solved off the grid (the observation-collar pass's workers): the solution and the gate's verdict, applied and
 * printed in the hole's turn. */
typedef struct CsInfillOut {
    float *pos, *tu, *tv;          /* [3 count] when kept */
    int gate;                      /* refused with the gate's line due */
    const char *reason;
    size_t failed_at;
    double failed_lo, failed_hi, failed_dot, failed_det;
} CsInfillOut;

static double cs_infill_coordinate(const CloseupGrid *g, const double *uv, size_t at, int d)
{
    return d < 3 ? g->pos[3*at+(size_t)d] : uv[2*at+(size_t)(d-3)];
}

static int cs_infill_known(const CloseupGrid *g, size_t at, int fill_collar)
{
    return g->state[at] == CLOSEUP_NODE_DATA ||
           (fill_collar && g->state[at] == CLOSEUP_NODE_FILL);
}

static void cs_infill_row(CsInfillSystem *s, const size_t *nodes, const double *weights, int count)
{
    if (s->nt < 0) return;
    int active = 0;
    double fixed[5] = {0};
    for (int k = 0; k < count; k++) {
        size_t at = nodes[k];
        if (cs_infill_idx(s, at) >= 0) { active++; continue; }
        if (!cs_infill_known(s->g,at,s->fill_collar) || s->g->component[at] != s->component) return;
        for (int d = 0; d < s->dimensions; d++) {
            double p = cs_infill_coordinate(s->g,s->uv,at,d);
            if (!isfinite(p)) return;
            fixed[d] += weights[k]*(p-s->origin[d]);
        }
    }
    if (!active) return;
    for (int k = 0; k < count; k++) {
        int a = cs_infill_idx(s, nodes[k]); if (a < 0) continue;
        for (int d = 0; d < s->dimensions; d++) s->rhs[(size_t)d*(size_t)s->n+(size_t)a] -= weights[k]*fixed[d];
        for (int j = 0; j <= k; j++) {
            int b = cs_infill_idx(s, nodes[j]); if (b < 0) continue;
            if (s->nt >= s->cap) { s->nt = -1; return; }
            s->rows[s->nt] = a > b ? a : b;
            s->cols[s->nt] = a > b ? b : a;
            s->values[s->nt++] = weights[k]*weights[j];
        }
    }
}

static double cs_infill_point(const CsInfillSystem *s, const double *solution, size_t at, int d)
{
    const int32_t a = cs_infill_idx(s, at);
    if (a < 0) return cs_infill_coordinate(s->g,s->uv,at,d);
    double value = solution[(size_t)d*(size_t)s->n+(size_t)a]+s->origin[d];
    /* The gates inspect the same float world points that will be exported. */
    return d < 3 ? (double)(float)value : value;
}

/* owner == NULL: index is set for the hole's nodes and cleared after, and a kept solution is written to the grid.
 * owner: index already holds every hole's numbering (owner[at] = its hole), and the solution and the gate's verdict
 * go to out (the grid is only read: several holes may solve at once). */
static int cs_infill_solve(CloseupGrid *g, const int32_t *nodes, size_t count, int32_t *index,
                            size_t x0, size_t y0, size_t x1, size_t y1, int32_t component,
                            const double normal[3], int fill_collar, int preserve_faces,
                            int provisional, double *uv, const int32_t *owner, int32_t region, CsInfillOut *out)
{
    size_t box = (x1-x0+3)*(y1-y0+3);
    if (count > CS_INFILL_MAX_NODES || box > 4*CS_INFILL_MAX_NODES) return 0;
    CsInfillSystem s = {0}; SparseFactor_T factor = NULL; double *solution = NULL;
    int kept = 0;
    const char *reason = "equations";
    size_t failed_at = 0;
    double failed_lo = 0, failed_hi = 0, failed_dot = 0, failed_det = 0;
    s.g = g; s.index = index; s.owner = owner; s.region = region; s.component = component; s.fill_collar = fill_collar;
    s.uv = uv; s.dimensions = uv ? 5 : 3;
    s.n = (int)count; s.cap = (int)(22*box);
    s.rows = (int *)malloc((size_t)s.cap*sizeof(int)); s.cols = (int *)malloc((size_t)s.cap*sizeof(int));
    s.values = (double *)malloc((size_t)s.cap*sizeof(double));
    s.rhs = (double *)calloc((size_t)s.dimensions*count,sizeof(double));
    solution = (double *)malloc((size_t)s.dimensions*count*sizeof(double));
    if (!s.rows || !s.cols || !s.values || !s.rhs || !solution) { kept = -1; goto done; }
    /* Local world origin avoids cancellation at large scan coordinates. */
    for (size_t j = y0-1; j <= y1+1; j++) for (size_t i = x0-1; i <= x1+1; i++) {
        size_t at = j*g->ni+i;
        if (cs_infill_known(g,at,fill_collar)) {
            for (int d = 0; d < s.dimensions; d++) s.origin[d] = cs_infill_coordinate(g,uv,at,d);
            goto origin_found;
        }
    }
origin_found:
    if (owner == NULL) for (size_t k = 0; k < count; k++) index[(size_t)nodes[k]] = (int32_t)k;
    for (size_t j = y0-1; j <= y1+1; j++) for (size_t i = x0-1; i <= x1+1; i++) {
        size_t at = j*g->ni+i;
        const double second[3] = {1,-2,1};
        const double mixed[4] = {.3535533905932738,-.3535533905932738,-.3535533905932738,.3535533905932738};
        /* ||X_uu||^2 + 2||X_uv||^2 + ||X_vv||^2, with a fixed two-row collar. */
        if (i > 0 && i+1 < g->ni) {
            size_t row[3] = {at-1,at,at+1}; cs_infill_row(&s,row,second,3);
        }
        if (j > 0 && j+1 < g->nj) {
            size_t row[3] = {at-g->ni,at,at+g->ni}; cs_infill_row(&s,row,second,3);
        }
        if (i > 0 && i+1 < g->ni && j > 0 && j+1 < g->nj) {
            size_t row[4] = {at-g->ni-1,at-g->ni+1,at+g->ni-1,at+g->ni+1};
            cs_infill_row(&s,row,mixed,4);
        }
        if (s.nt < 0) goto done;
    }
    reason = "factor";
    if (!s.nt || Sparse_factor_spd(s.n,s.nt,s.rows,s.cols,s.values,&factor) ||
        Sparse_factor_solve_multi(factor,s.rhs,solution,s.dimensions)) goto done;
    /* Test every triangle touching this hole, including its diagonal and
     * boundary edges. A folded or stretched solution is refused as a whole. */
    for (size_t j = y0-1; j <= y1; j++) for (size_t i = x0-1; i <= x1; i++) {
        size_t a = j*g->ni+i, corners[4] = {a,a+1,a+g->ni,a+g->ni+1};
        int tri[2][3] = {{0,1,3},{0,3,2}};
        if (preserve_faces) {
            int unknown = 0, missing = -1;
            for (int k = 0; k < 4; k++) if (cs_infill_idx(&s, corners[k]) >= 0) { unknown++; missing = k; }
            /* A retained three-corner face owns its old diagonal. Inserting
             * the fourth corner adds its complement, never retriangulates or
             * overlaps that face. The exporter must use the same topology. */
            if (unknown == 1 && (missing == 0 || missing == 3)) {
                tri[0][0] = 0; tri[0][1] = 1; tri[0][2] = 2;
                tri[1][0] = 1; tri[1][1] = 3; tri[1][2] = 2;
            }
        }
        for (int t = 0; t < 2; t++) {
            size_t v[3]; int active = 0, valid = 1;
            double e[3], f[3], cross[3], ee = 0, ff = 0, ef = 0, dot = 0;
            for (int k = 0; k < 3; k++) {
                v[k] = corners[tri[t][k]]; active += cs_infill_idx(&s, v[k]) >= 0;
                if (cs_infill_idx(&s, v[k]) < 0 && (!cs_infill_known(g,v[k],fill_collar) || g->component[v[k]] != component)) valid = 0;
            }
            if (!active || !valid) continue;
            double au = g->h*((double)(v[1]%g->ni)-(double)(v[0]%g->ni));
            double av = g->h*((double)(v[1]/g->ni)-(double)(v[0]/g->ni));
            double bu = g->h*((double)(v[2]%g->ni)-(double)(v[0]%g->ni));
            double bv = g->h*((double)(v[2]/g->ni)-(double)(v[0]/g->ni));
            if (uv) {
                au = cs_infill_point(&s,solution,v[1],3)-cs_infill_point(&s,solution,v[0],3);
                av = cs_infill_point(&s,solution,v[1],4)-cs_infill_point(&s,solution,v[0],4);
                bu = cs_infill_point(&s,solution,v[2],3)-cs_infill_point(&s,solution,v[0],3);
                bv = cs_infill_point(&s,solution,v[2],4)-cs_infill_point(&s,solution,v[0],4);
            }
            double det = au*bv-av*bu;
            if (!(det > 0) || !isfinite(det)) {
                reason = "uv-orientation"; failed_at = a; failed_det = det; goto done;
            }
            for (int d = 0; d < 3; d++) {
                double p = cs_infill_point(&s,solution,v[0],d);
                double edge = cs_infill_point(&s,solution,v[1],d)-p;
                double across = cs_infill_point(&s,solution,v[2],d)-p;
                e[d] = (edge*bv-across*av)/det;
                f[d] = (across*au-edge*bu)/det;
                ee += e[d]*e[d]; ff += f[d]*f[d]; ef += e[d]*f[d];
            }
            for (int d = 0; d < 3; d++) {
                cross[d] = e[(d+1)%3]*f[(d+2)%3]-e[(d+2)%3]*f[(d+1)%3];
                dot += cross[d]*normal[d];
            }
            double disc = hypot(ee-ff,2*ef), hi = .5*(ee+ff+disc);
            double area2 = cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2];
            double lo = hi > 0 ? area2/hi : 0;
            if (!isfinite(hi) || !isfinite(lo) || !(lo > 0) || !(dot > 0) ||
                (!provisional && (lo < .25 || hi > 4))) {
                reason = !isfinite(hi) || !isfinite(lo) ? "nonfinite" : !(lo > 0) ? "world-degenerate" :
                         !(dot > 0) ? "orientation" : lo < .25 ? "compression" : "stretch";
                failed_at = a; failed_lo = lo; failed_hi = hi; failed_dot = dot; failed_det = det;
                goto done;
            }
        }
    }
    if (out != NULL) {          /* the same values, for the hole's turn (cs_fill_holes_observed) */
        out->pos = (float *)malloc(3*count*sizeof(float)); out->tu = (float *)malloc(3*count*sizeof(float));
        out->tv = (float *)malloc(3*count*sizeof(float));
        if (!out->pos || !out->tu || !out->tv) { kept = -1; goto done; }
        for (size_t k = 0; k < count; k++) {
            size_t at = (size_t)nodes[k];
            for (int d = 0; d < 3; d++) {
                out->pos[3*k+(size_t)d] = (float)cs_infill_point(&s,solution,at,d);
                out->tu[3*k+(size_t)d] = (float)((cs_infill_point(&s,solution,at+1,d)-
                    cs_infill_point(&s,solution,at-1,d))/(2*g->h));
                out->tv[3*k+(size_t)d] = (float)((cs_infill_point(&s,solution,at+g->ni,d)-
                    cs_infill_point(&s,solution,at-g->ni,d))/(2*g->h));
            }
        }
        kept = 1;
        goto done;
    }
    for (size_t k = 0; k < count; k++) {
        size_t at = (size_t)nodes[k];
        for (int d = 0; d < 3; d++) {
            g->pos[3*at+(size_t)d] = (float)cs_infill_point(&s,solution,at,d);
            g->tu[3*at+(size_t)d] = (float)((cs_infill_point(&s,solution,at+1,d)-
                cs_infill_point(&s,solution,at-1,d))/(2*g->h));
            g->tv[3*at+(size_t)d] = (float)((cs_infill_point(&s,solution,at+g->ni,d)-
                cs_infill_point(&s,solution,at-g->ni,d))/(2*g->h));
        }
        if (uv) for (int d = 0; d < 2; d++) uv[2*at+(size_t)d] = cs_infill_point(&s,solution,at,d+3);
        g->state[at] = CLOSEUP_NODE_FILL; g->chart[at] = -1; g->component[at] = component;
        g->support[at] = 0; cs_set_normal(g,at);
    }
    kept = 1;
done:
    if (!kept && out != NULL) {
        out->gate = 1; out->reason = reason; out->failed_at = failed_at;
        out->failed_lo = failed_lo; out->failed_hi = failed_hi; out->failed_dot = failed_dot; out->failed_det = failed_det;
    } else if (!kept) fprintf(stderr,"    Cartesian infill gate: %s box=%zu,%zu..%zu,%zu at=%zu,%zu "
                      "squared-scales=%.8g,%.8g orientation=%.8g uv-det=%.8g\n",reason,x0,y0,x1,y1,
                      failed_at%g->ni,failed_at/g->ni,failed_lo,failed_hi,failed_dot,failed_det);
    if (owner == NULL) for (size_t k = 0; k < count; k++) index[(size_t)nodes[k]] = -1;
    Sparse_factor_free(&factor); free(s.rows); free(s.cols); free(s.values); free(s.rhs); free(solution);
    return kept;
}

static int cs_fill_holes_pass(CloseupGrid *g, double max_distance_px, CloseupFitStats *st,
                               int fill_collar, int require_fill, int preserve_faces,
                               int provisional, int32_t *patch_ids,
                               const float *reference_normal, double *uv)
{
    if (!g || !st || !(max_distance_px > 0) || !isfinite(max_distance_px) || !(g->h > 0) ||
        !g->ni || !g->nj || g->ni > INT32_MAX/g->nj) return -1;
    size_t nn = g->ni*g->nj;
    uint8_t *seen = (uint8_t *)calloc(nn,1);
    int32_t *queue = (int32_t *)malloc(nn*sizeof(int32_t)), *index = (int32_t *)malloc(nn*sizeof(int32_t));
    int rc = -1;
    double t_pass = ves_clock_sec(), t_solve = 0;
    size_t n_solve = 0, n_solve_nodes = 0, n_regions = 0;
    if (!seen || !queue || !index) goto done;
    for (size_t k = 0; k < nn; k++) index[k] = -1;
    for (size_t start = 0; start < nn; start++) {
        if (seen[start] || cs_infill_known(g,start,fill_collar)) continue;
        size_t count = 1, x0 = g->ni, x1 = 0, y0 = g->nj, y1 = 0, added = 0;
        int exterior = 0, mixed = 0, beyond_limit = 0, has_fill = 0;
        int32_t component = -1; double normal[3] = {0};
        seen[start] = 1; queue[0] = (int32_t)start;
        for (size_t head = 0; head < count; head++) {
            size_t at = (size_t)queue[head], x = at%g->ni, y = at/g->ni;
            if (!x || !y || x+1 == g->ni || y+1 == g->nj) exterior = 1;
            if (x < x0) x0 = x; if (x > x1) x1 = x;
            if (y < y0) y0 = y; if (y > y1) y1 = y;
            beyond_limit |= !(g->dist[at] <= max_distance_px);
            added += g->state[at] != CLOSEUP_NODE_FILL;
            for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
                if ((!dx && !dy) || (dx < 0 && !x) || (dy < 0 && !y) ||
                    (dx > 0 && x+1 == g->ni) || (dy > 0 && y+1 == g->nj)) continue;
                size_t other = (size_t)((int64_t)at+(int64_t)dy*(int64_t)g->ni+dx);
                if (cs_infill_known(g,other,fill_collar)) {
                    has_fill |= g->state[other] == CLOSEUP_NODE_FILL;
                    int32_t c = g->component[other];
                    if (c < 0) mixed = 1;
                    if (component < 0) component = c; else if (component != c) mixed = 1;
                    for (int d = 0; d < 3; d++) normal[d] += reference_normal[3*other+(size_t)d];
                } else if (!seen[other]) { seen[other] = 1; queue[count++] = (int32_t)other; }
            }
        }
        n_regions++;
        if (exterior || (require_fill && !has_fill)) continue;
        st->variational_holes++;
        int result = 0;
        if (!mixed && !beyond_limit && component >= 0 && count <= CS_INFILL_MAX_NODES) {
            double ts = ves_clock_sec();
            result = cs_infill_solve(g,queue,count,index,x0,y0,x1,y1,component,normal,
                                    fill_collar,preserve_faces,provisional,uv,NULL,0,NULL);
            t_solve += ves_clock_sec() - ts; n_solve++; n_solve_nodes += count;
        }
        if (result < 0) goto done;
        if (result) {
            if (patch_ids) for (size_t k = 0; k < count; k++) patch_ids[(size_t)queue[k]] = (int32_t)st->variational_filled;
            st->variational_filled++; st->variational_nodes += count; st->fill_nodes += added;
        }
        else {
            st->variational_refused++;
            fprintf(stderr,"  Cartesian infill refused: collar=%s nodes=%zu box=%zu,%zu..%zu,%zu mixed=%d beyond-reach=%d\n",
                    fill_collar?"surface":"observations",count,x0,y0,x1,y1,mixed,beyond_limit);
        }
    }
    rc = 0;
done:
    fprintf(stderr, "  [holes-timing] collar=%d regions %zu solves %zu (%zu nodes) solve %.2fs pass %.2fs\n",
            fill_collar, n_regions, n_solve, n_solve_nodes, t_solve, ves_clock_sec() - t_pass);
    free(seen); free(queue); free(index); return rc;
}

/* cs_fill_holes_pass(g, max, st, 0, 0, 0, 0, NULL, reference_normal, NULL) with its holes solved on all threads.
 * With the observations alone as collar no hole reads another: a hole's rows and gates read only DATA nodes and its
 * own, and a filled hole is FILL, which this pass does not count as known.  So the holes are found in the scan's order,
 * solved at once off the grid, then applied and logged in that order: the same grid, counts and log. */
typedef struct CsHole {
    size_t first, count, x0, y0, x1, y1, added;
    int exterior, mixed, beyond, solvable;
    int32_t component;
    double normal[3];
} CsHole;

static int cs_fill_holes_observed(CloseupGrid *g, double max_distance_px, CloseupFitStats *st,
                                  const float *reference_normal)
{
    if (!g || !st || !(max_distance_px > 0) || !isfinite(max_distance_px) || !(g->h > 0) ||
        !g->ni || !g->nj || g->ni > INT32_MAX/g->nj) return -1;
    size_t nn = g->ni*g->nj, used = 0, nh = 0, caph = 1024;
    uint8_t *seen = (uint8_t *)calloc(nn,1);
    int32_t *nodes = (int32_t *)malloc(nn*sizeof(int32_t)), *index = (int32_t *)malloc(nn*sizeof(int32_t));
    int32_t *owner = (int32_t *)malloc(nn*sizeof(int32_t));
    CsHole *holes = (CsHole *)malloc(caph*sizeof *holes);
    CsInfillOut *res = NULL;
    int *kept = NULL, rc = -1, bad = 0;
    long r = 0;
    double t_pass = ves_clock_sec();
    size_t n_solve = 0, n_solve_nodes = 0;
    if (!seen || !nodes || !index || !owner || !holes) goto done;
    for (size_t k = 0; k < nn; k++) { index[k] = -1; owner[k] = -1; }
    for (size_t start = 0; start < nn; start++) {   /* the holes, as the serial scan meets them */
        if (seen[start] || cs_infill_known(g,start,0)) continue;
        CsHole h;
        memset(&h, 0, sizeof h);
        h.first = used; h.count = 1; h.x0 = g->ni; h.y0 = g->nj; h.component = -1;
        seen[start] = 1; nodes[used] = (int32_t)start;
        for (size_t head = 0; head < h.count; head++) {
            size_t at = (size_t)nodes[used+head], x = at%g->ni, y = at/g->ni;
            if (!x || !y || x+1 == g->ni || y+1 == g->nj) h.exterior = 1;
            if (x < h.x0) h.x0 = x; if (x > h.x1) h.x1 = x;
            if (y < h.y0) h.y0 = y; if (y > h.y1) h.y1 = y;
            h.beyond |= !(g->dist[at] <= max_distance_px);
            h.added += g->state[at] != CLOSEUP_NODE_FILL;
            for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
                if ((!dx && !dy) || (dx < 0 && !x) || (dy < 0 && !y) ||
                    (dx > 0 && x+1 == g->ni) || (dy > 0 && y+1 == g->nj)) continue;
                size_t other = (size_t)((int64_t)at+(int64_t)dy*(int64_t)g->ni+dx);
                if (cs_infill_known(g,other,0)) {
                    int32_t c = g->component[other];
                    if (c < 0) h.mixed = 1;
                    if (h.component < 0) h.component = c; else if (h.component != c) h.mixed = 1;
                    for (int d = 0; d < 3; d++) h.normal[d] += reference_normal[3*other+(size_t)d];
                } else if (!seen[other]) { seen[other] = 1; nodes[used+h.count++] = (int32_t)other; }
            }
        }
        used += h.count;
        h.solvable = !h.exterior && !h.mixed && !h.beyond && h.component >= 0 && h.count <= CS_INFILL_MAX_NODES;
        if (nh == caph) {
            CsHole *grown = (CsHole *)realloc(holes, 2*caph*sizeof *holes);
            if (!grown) goto done;
            holes = grown; caph *= 2;
        }
        holes[nh++] = h;
    }
    res = (CsInfillOut *)calloc(nh ? nh : 1, sizeof *res);
    kept = (int *)calloc(nh ? nh : 1, sizeof *kept);
    if (!res || !kept) goto done;
    for (size_t h = 0; h < nh; h++) {
        if (!holes[h].solvable) continue;
        for (size_t k = 0; k < holes[h].count; k++) {
            index[nodes[holes[h].first+k]] = (int32_t)k;
            owner[nodes[holes[h].first+k]] = (int32_t)h;
        }
        n_solve++; n_solve_nodes += holes[h].count;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(| : bad)
#endif
    for (r = 0; r < (long)nh; r++) {
        const CsHole *h = &holes[r];
        if (!h->solvable) continue;
        kept[r] = cs_infill_solve(g,nodes+h->first,h->count,index,h->x0,h->y0,h->x1,h->y1,h->component,h->normal,
                                  0,0,0,NULL,owner,(int32_t)r,&res[r]);
        if (kept[r] < 0) bad |= 1;
    }
    if (bad) goto done;
    for (size_t h = 0; h < nh; h++) {   /* their turns: the serial pass's writes, counts and log */
        const CsHole *o = &holes[h];
        if (o->exterior) continue;
        st->variational_holes++;
        if (o->solvable && !kept[h] && res[h].gate)
            fprintf(stderr,"    Cartesian infill gate: %s box=%zu,%zu..%zu,%zu at=%zu,%zu "
                    "squared-scales=%.8g,%.8g orientation=%.8g uv-det=%.8g\n",res[h].reason,o->x0,o->y0,o->x1,o->y1,
                    res[h].failed_at%g->ni,res[h].failed_at/g->ni,res[h].failed_lo,res[h].failed_hi,
                    res[h].failed_dot,res[h].failed_det);
        if (o->solvable && kept[h]) {
            for (size_t k = 0; k < o->count; k++) {
                size_t at = (size_t)nodes[o->first+k];
                memcpy(&g->pos[3*at], &res[h].pos[3*k], 3*sizeof(float));
                memcpy(&g->tu[3*at], &res[h].tu[3*k], 3*sizeof(float));
                memcpy(&g->tv[3*at], &res[h].tv[3*k], 3*sizeof(float));
                g->state[at] = CLOSEUP_NODE_FILL; g->chart[at] = -1; g->component[at] = o->component;
                g->support[at] = 0; cs_set_normal(g,at);
            }
            st->variational_filled++; st->variational_nodes += o->count; st->fill_nodes += o->added;
        } else {
            st->variational_refused++;
            fprintf(stderr,"  Cartesian infill refused: collar=%s nodes=%zu box=%zu,%zu..%zu,%zu mixed=%d beyond-reach=%d\n",
                    "observations",o->count,o->x0,o->y0,o->x1,o->y1,o->mixed,o->beyond);
        }
    }
    rc = 0;
done:
    fprintf(stderr, "  [holes-timing] collar=0 (parallel) regions %zu solves %zu (%zu nodes) pass %.2fs\n",
            nh, n_solve, n_solve_nodes, ves_clock_sec() - t_pass);
    if (res) for (size_t h = 0; h < nh; h++) { free(res[h].pos); free(res[h].tu); free(res[h].tv); }
    free(seen); free(nodes); free(index); free(owner); free(holes); free(res); free(kept);
    return rc;
}

int CloseupSurface_fill_holes(CloseupGrid *g, double max_distance_px, CloseupFitStats *st)
{
    if (!g) return -1;
    int rc=cs_fill_holes_observed(g,max_distance_px,st,g->nrm);
    return rc?rc:cs_fill_holes_pass(g,max_distance_px,st,1,1,0,0,NULL,g->nrm,NULL);
}

/* Importers retain their real face-component labels. Mere grid adjacency must
 * never merge two wraps. The geometric orientation reference is separate from
 * the stored shading normals, which may have been turned inward. */
static void cs_infill_reference_triangle(const CloseupGrid *g, size_t a, size_t b, size_t c,
                                          float *reference)
{
    if (!g->state[a] || !g->state[b] || !g->state[c] || g->component[a] < 0 ||
        g->component[a] != g->component[b] || g->component[a] != g->component[c]) return;
    double e[3], f[3], n[3], length = 0;
    for (int d = 0; d < 3; d++) {
        e[d] = (double)g->pos[3*b+(size_t)d]-g->pos[3*a+(size_t)d];
        f[d] = (double)g->pos[3*c+(size_t)d]-g->pos[3*a+(size_t)d];
    }
    for (int d = 0; d < 3; d++) {
        n[d] = e[(d+1)%3]*f[(d+2)%3]-e[(d+2)%3]*f[(d+1)%3];
        length += n[d]*n[d];
    }
    if (!(length > 0) || !isfinite(length)) return;
    /* Unit triangle votes keep one stretched retained face from deciding the
     * sign of a whole hole. This reference never changes retained geometry. */
    length = sqrt(length);
    for (int d = 0; d < 3; d++) {
        float value = (float)(n[d]/length);
        reference[3*a+(size_t)d] += value;
        reference[3*b+(size_t)d] += value;
        reference[3*c+(size_t)d] += value;
    }
}

static int cs_complete_fixed(CloseupGrid *g, double *uv_nodes, double max_distance_px,
                               CloseupFitStats *st, int provisional, int32_t *patch_ids)
{
    if (!g || !st || !g->ni || !g->nj || g->ni > INT32_MAX/g->nj || !(g->h > 0) ||
        !(max_distance_px > 0) || !isfinite(max_distance_px)) return -1;
    size_t nn = g->ni*g->nj;
    uint8_t *seed = (uint8_t *)calloc(nn,1), *was_valid = (uint8_t *)calloc(nn,1);
    float *d2 = (float *)malloc(nn*sizeof(float)), *reference = (float *)calloc(3*nn,sizeof(float));
    int32_t *nearest = (int32_t *)malloc(nn*sizeof(int32_t));
    int rc = -1;
    if (!seed || !was_valid || !d2 || !reference || !nearest) goto done;
    for (size_t k = 0; k < nn; k++) {
        if (!g->state[k]) continue;
        if ((g->state[k] != CLOSEUP_NODE_DATA && g->state[k] != CLOSEUP_NODE_FILL) || g->component[k] < 0)
            goto done;
        for (int d = 0; d < 3; d++) if (!isfinite(g->pos[3*k+(size_t)d])) goto done;
        if (uv_nodes && (!isfinite(uv_nodes[2*k]) || !isfinite(uv_nodes[2*k+1]))) goto done;
        was_valid[k] = 1; seed[k] = g->state[k] == CLOSEUP_NODE_DATA;
    }
    if (cs_edt2d(seed,g->ni,g->nj,d2,nearest)) goto done;
    for (size_t k = 0; k < nn; k++) if (!was_valid[k])
        g->dist[k] = d2[k] < CS_INF ? (float)(g->h*sqrt((double)d2[k])) : INFINITY;
    for (size_t j = 0; j+1 < g->nj; j++) for (size_t i = 0; i+1 < g->ni; i++) {
        size_t a = j*g->ni+i, b = a+1, c = a+g->ni, d = c+1;
        if (g->state[a] && g->state[b] && g->state[c] && g->state[d]) {
            cs_infill_reference_triangle(g,a,b,d,reference); cs_infill_reference_triangle(g,a,d,c,reference);
        } else if (g->state[a] && g->state[b] && g->state[c]) cs_infill_reference_triangle(g,a,b,c,reference);
        else if (g->state[b] && g->state[d] && g->state[c]) cs_infill_reference_triangle(g,b,d,c,reference);
        else if (g->state[a] && g->state[d] && g->state[c]) cs_infill_reference_triangle(g,a,d,c,reference);
        else if (g->state[a] && g->state[b] && g->state[d]) cs_infill_reference_triangle(g,a,b,d,reference);
    }
    for (size_t k = 0; k < nn; k++) {
        double length = hypot(hypot(reference[3*k],reference[3*k+1]),reference[3*k+2]);
        if (length > 0) for (int d = 0; d < 3; d++) reference[3*k+(size_t)d] = (float)(reference[3*k+(size_t)d]/length);
    }
    rc = cs_fill_holes_pass(g,max_distance_px,st,1,0,1,provisional,patch_ids,reference,uv_nodes);
    /* Newly synthesized shading normals inherit the saved surface's inward
     * sign; the orientation gate itself uses UV-oriented geometry. */
    for (size_t k = 0; k < nn; k++) {
        if (was_valid[k] || !g->state[k] || nearest[k] < 0) continue;
        size_t q = (size_t)nearest[k]; double dot = 0;
        for (int d = 0; d < 3; d++) dot += (double)g->nrm[3*q+(size_t)d]*reference[3*q+(size_t)d];
        if (dot < 0) for (int d = 0; d < 3; d++) g->nrm[3*k+(size_t)d] = -g->nrm[3*k+(size_t)d];
    }
done:
    free(seed); free(was_valid); free(d2); free(reference); free(nearest);
    return rc;
}

int CloseupSurface_complete_fixed(CloseupGrid *g, double max_distance_px, CloseupFitStats *st)
{
    return CloseupSurface_complete_fixed_uv(g,NULL,max_distance_px,st);
}

int CloseupSurface_complete_fixed_uv(CloseupGrid *g, double *uv_nodes, double max_distance_px, CloseupFitStats *st)
{
    return cs_complete_fixed(g,uv_nodes,max_distance_px,st,0,NULL);
}

#include "closeup_provisional.inc"
#include "closeup_rebuild.inc"

static int cs_selftest_fixed_infill(void)
{
    int fails = 0;
    for (int test = 0; test < 7; test++) {
        CloseupGrid g = {0}; CloseupFitStats st = {0};
        if (CloseupGrid_alloc(&g,-40,40,-40,40,2)) return 1;
        size_t nn = g.ni*g.nj, hole = 0, filled = 0;
        float *saved = (float *)malloc(14*nn*sizeof(float));
        int32_t *labels = (int32_t *)malloc(2*nn*sizeof(int32_t));
        uint8_t *state = (uint8_t *)malloc(nn);
        double *uv = (double *)malloc(4*nn*sizeof(double));
        if (!saved || !labels || !state || !uv) {
            free(saved); free(labels); free(state); free(uv); CloseupGrid_free(&g); return 1;
        }
        for (size_t k = 0; k < nn; k++) {
            double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
            double scale = test == 3 || test == 5 ? .05 : 1;
            g.pos[3*k] = (float)(12000+scale*u); g.pos[3*k+1] = (float)(8000+scale*v);
            g.pos[3*k+2] = (float)(4000+.02*sin(u/10)*sin(v/7));
            g.tu[3*k] = (float)scale; g.tv[3*k+1] = (float)scale;
            g.tu[3*k+2] = (float)(.002*cos(u/10)*sin(v/7));
            g.tv[3*k+2] = (float)(.02/7*sin(u/10)*cos(v/7));
            cs_set_normal(&g,k);
            if (test == 1) for (int d = 0; d < 3; d++) g.nrm[3*k+(size_t)d] = -g.nrm[3*k+(size_t)d];
            g.state[k] = test != 0 && fabs(u)<32 && fabs(v)<32 ? CLOSEUP_NODE_FILL : CLOSEUP_NODE_DATA;
            g.component[k] = test == 2 && u>0 ? 2 : 1;
            g.chart[k] = g.state[k] == CLOSEUP_NODE_DATA ? 17 : -1;
            g.dist[k] = g.state[k] == CLOSEUP_NODE_DATA ? 0 : 6;
            g.support[k] = g.state[k] == CLOSEUP_NODE_DATA ? 12 : 0;
            uv[2*k] = 1e6+(test == 6 ? -1 : scale)*u;
            uv[2*k+1] = 2e6+scale*v;
            if ((fabs(u)<8 && fabs(v)<10) || (u>30 && fabs(v)<6)) {
                g.state[k] = CLOSEUP_NODE_NONE; g.component[k] = -1; g.chart[k] = -1;
                for (int d = 0; d < 3; d++) g.pos[3*k+(size_t)d] = NAN;
                uv[2*k] = uv[2*k+1] = NAN;
                if (fabs(u)<8) hole++;
            }
        }
        float *arrays[6] = {g.pos,g.nrm,g.tu,g.tv,g.dist,g.support};
        const size_t widths[6] = {3,3,3,3,1,1}; size_t offset = 0;
        for (int a = 0; a < 6; a++) {
            memcpy(saved+offset,arrays[a],widths[a]*nn*sizeof(float)); offset += widths[a]*nn;
        }
        memcpy(labels,g.chart,nn*sizeof(int32_t)); memcpy(labels+nn,g.component,nn*sizeof(int32_t));
        memcpy(state,g.state,nn);
        memcpy(uv+2*nn,uv,2*nn*sizeof(double));
        if (test >= 5 ? CloseupSurface_complete_fixed_uv(&g,uv,64,&st) :
                       CloseupSurface_complete_fixed(&g,test == 4 ? 2 : 64,&st)) fails++;
        for (size_t k = 0; k < nn; k++) {
            if (state[k]) {
                if (g.state[k] != state[k] || g.chart[k] != labels[k] || g.component[k] != labels[nn+k]) fails++;
                if (memcmp(uv+2*k,uv+2*nn+2*k,2*sizeof(double))) fails++;
                offset = 0;
                for (int a = 0; a < 6; a++) {
                    if (memcmp(arrays[a]+widths[a]*k,saved+offset+widths[a]*k,widths[a]*sizeof(float))) fails++;
                    offset += widths[a]*nn;
                }
            } else {
                double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
                int should_fill = (test<2 || test == 5) && fabs(u)<8 && fabs(v)<10;
                if (should_fill) {
                    filled += g.state[k] == CLOSEUP_NODE_FILL;
                    if (g.state[k] != CLOSEUP_NODE_FILL || (test == 1 ? g.nrm[3*k+2] >= 0 : g.nrm[3*k+2] <= 0)) fails++;
                    double scale = test == 5 ? .05 : 1;
                    if (fabs((double)g.pos[3*k]-(12000+scale*u))>.01 ||
                        fabs((double)g.pos[3*k+1]-(8000+scale*v))>.01) fails++;
                    if (test == 5 && (fabs(uv[2*k]-(1e6+.05*u))>1e-7 || fabs(uv[2*k+1]-(2e6+.05*v))>1e-7)) fails++;
                } else {
                    if (g.state[k]) fails++;
                    if (memcmp(uv+2*k,uv+2*nn+2*k,2*sizeof(double))) fails++;
                }
            }
        }
        if (test<2 || test == 5 ? filled != hole || st.variational_filled != 1 || st.variational_nodes != hole :
                     st.variational_filled != 0 || st.variational_refused != 1) fails++;
        fprintf(stderr,"  Fixed Cartesian infill case %d: %zu filled, %zu refused; retained arrays exact\n",
                test,st.variational_filled,st.variational_refused);
        free(saved); free(labels); free(state); free(uv); CloseupGrid_free(&g);
    }
    fprintf(stderr,"  Fixed Cartesian infill: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails != 0;
}

/* The parallel observation-collar pass against the serial one on copies of a sheet with ~300 holes: small and large,
 * kept and gated, mixed components, beyond reach, at the border; every array and counter must agree. */
static int cs_selftest_infill_parallel(void)
{
    CloseupGrid ga = {0}, gb = {0};
    CloseupFitStats sa = {0}, sb = {0};
    uint32_t rng = 0x3c6ef372u;
    size_t nn = 0, k = 0, bad = 0;
    int fails = 0, h = 0;
    if (CloseupGrid_alloc(&ga,0,600,0,400,2) || CloseupGrid_alloc(&gb,0,600,0,400,2)) {
        CloseupGrid_free(&ga); CloseupGrid_free(&gb); return 1;
    }
    nn = ga.ni*ga.nj;
    for (k = 0; k < nn; k++) {
        double u = ga.u0+ga.h*(double)(k%ga.ni), v = ga.v0+ga.h*(double)(k/ga.ni);
        ga.pos[3*k] = (float)(9000+u); ga.pos[3*k+1] = (float)(7000+v);
        ga.pos[3*k+2] = (float)(3000+.02*u*v/40+3*sin(u/37));
        ga.tu[3*k] = 1; ga.tv[3*k+1] = 1;
        ga.tu[3*k+2] = (float)(.02*v/40+3*cos(u/37)/37); ga.tv[3*k+2] = (float)(.02*u/40);
        cs_set_normal(&ga,k);
        ga.state[k] = CLOSEUP_NODE_DATA; ga.component[k] = u > 400 ? 2 : 1; ga.chart[k] = 5; ga.dist[k] = 0;
        ga.support[k] = 3;
    }
    for (h = 0; h < 300; h++) {        /* rectangles of 1..12 x 1..9 nodes; some stretched by a bump */
        size_t x0 = 0, y0 = 0, w = 0, hh = 0, x = 0, y = 0;
        rng = rng*1664525u+1013904223u; x0 = (rng >> 8) % ga.ni;
        rng = rng*1664525u+1013904223u; y0 = (rng >> 8) % ga.nj;
        rng = rng*1664525u+1013904223u; w = 1 + (rng >> 28) % 12; hh = 1 + (rng >> 24) % 9;
        for (y = y0; y < y0+hh && y < ga.nj; y++) for (x = x0; x < x0+w && x < ga.ni; x++) {
            size_t at = y*ga.ni+x;
            ga.state[at] = (rng >> 20) % 5 == 0 ? CLOSEUP_NODE_FILL : CLOSEUP_NODE_NONE;
            ga.component[at] = ga.state[at] == CLOSEUP_NODE_FILL ? ga.component[at] : -1;
            ga.dist[at] = (rng >> 16) % 23 == 0 ? 1000.0f : 4.0f;
            if (ga.state[at] == CLOSEUP_NODE_NONE) ga.pos[3*at] = ga.pos[3*at+1] = ga.pos[3*at+2] = NAN;
            if ((rng >> 12) % 7 == 0) ga.pos[3*at+2] += 30.0f;
        }
    }
    gb.ni = ga.ni; gb.nj = ga.nj; gb.u0 = ga.u0; gb.v0 = ga.v0; gb.h = ga.h;
    memcpy(gb.pos,ga.pos,3*nn*sizeof(float)); memcpy(gb.tu,ga.tu,3*nn*sizeof(float)); memcpy(gb.tv,ga.tv,3*nn*sizeof(float));
    memcpy(gb.nrm,ga.nrm,3*nn*sizeof(float)); memcpy(gb.dist,ga.dist,nn*sizeof(float));
    memcpy(gb.support,ga.support,nn*sizeof(float)); memcpy(gb.chart,ga.chart,nn*sizeof(int32_t));
    memcpy(gb.component,ga.component,nn*sizeof(int32_t)); memcpy(gb.state,ga.state,nn);
    if (cs_fill_holes_pass(&ga,64,&sa,0,0,0,0,NULL,ga.nrm,NULL) || cs_fill_holes_observed(&gb,64,&sb,gb.nrm)) fails++;
    for (k = 0; k < nn; k++)
        if (ga.state[k] != gb.state[k] || ga.chart[k] != gb.chart[k] || ga.component[k] != gb.component[k] ||
            memcmp(&ga.pos[3*k],&gb.pos[3*k],3*sizeof(float)) || memcmp(&ga.tu[3*k],&gb.tu[3*k],3*sizeof(float)) ||
            memcmp(&ga.tv[3*k],&gb.tv[3*k],3*sizeof(float)) || memcmp(&ga.nrm[3*k],&gb.nrm[3*k],3*sizeof(float)) ||
            memcmp(&ga.support[k],&gb.support[k],sizeof(float))) bad++;
    if (sa.variational_holes != sb.variational_holes || sa.variational_filled != sb.variational_filled ||
        sa.variational_refused != sb.variational_refused || sa.variational_nodes != sb.variational_nodes ||
        sa.fill_nodes != sb.fill_nodes) bad++;
    if (bad || sa.variational_filled == 0 || sa.variational_refused == 0) fails++;
    fprintf(stderr,"  parallel observation-collar holes vs the serial pass: %zu holes, %zu filled, %zu refused, %zu "
            "differences: %s\n",sa.variational_holes,sa.variational_filled,sa.variational_refused,bad,fails ? "FAIL" : "ok");
    CloseupGrid_free(&ga); CloseupGrid_free(&gb);
    return fails;
}

static int cs_selftest_infill(void)
{
    int fails = 0;
    for (int test = 0; test < 7; test++) {
        CloseupGrid g = {0}; CloseupFitStats st = {0};
        if (CloseupGrid_alloc(&g,-64,64,-64,64,2)) return 1;
        size_t nn = g.ni*g.nj, hole = 0; double worst = 0;
        float *truth = (float *)malloc(3*nn*sizeof(float));
        if (!truth) { CloseupGrid_free(&g); return 1; }
        for (size_t k = 0; k < nn; k++) {
            double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
            double z = test == 1 ? 70*cos(u/70) : .003*u*u-.002*v*v;
            double x = test == 1 ? 70*sin(u/70) : u, scale = test == 3 || test == 6 ? .05 : 1;
            /* A proper world rotation plus a large scan translation. */
            g.pos[3*k] = (float)(12000+scale*(.8*z+.6*x));
            g.pos[3*k+1] = (float)(8000+scale*v);
            g.pos[3*k+2] = (float)(4000+scale*(-.6*z+.8*x));
            double zu = test == 1 ? -sin(u/70) : .006*u;
            double xu = test == 1 ? cos(u/70) : 1, zv = test == 1 ? 0 : -.004*v;
            g.tu[3*k] = (float)(scale*(.8*zu+.6*xu)); g.tu[3*k+1] = 0;
            g.tu[3*k+2] = (float)(scale*(-.6*zu+.8*xu));
            g.tv[3*k] = (float)(scale*.8*zv); g.tv[3*k+1] = (float)scale;
            g.tv[3*k+2] = (float)(-scale*.6*zv); cs_set_normal(&g,k);
            g.state[k] = test>=4 && fabs(u)<62 && fabs(v)<62 ? CLOSEUP_NODE_FILL : CLOSEUP_NODE_DATA;
            g.component[k] = (test == 2 || test == 5) && u > 0 ? 2 : 1;
            g.dist[k] = 0;
        }
        memcpy(truth,g.pos,3*nn*sizeof(float));
        for (size_t k = 0; k < nn; k++) {
            double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
            if ((fabs(u) < 14 && fabs(v) < 18) || (u > 50 && fabs(v) < 10)) {
                g.state[k] = CLOSEUP_NODE_NONE; g.component[k] = -1; g.dist[k] = 16;
                g.pos[3*k] = g.pos[3*k+1] = g.pos[3*k+2] = NAN;
                if (fabs(u) < 14) hole++;
            }
        }
        if (CloseupSurface_fill_holes(&g,64,&st)) fails++;
        for (size_t k = 0; k < nn; k++) {
            double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
            if (fabs(u) < 14 && fabs(v) < 18) {
                if (test < 2 || test == 4) {
                    if (g.state[k] != CLOSEUP_NODE_FILL) { fails++; continue; }
                    double e = 0;
                    for (int d = 0; d < 3; d++) {
                        double diff = g.pos[3*k+(size_t)d]-truth[3*k+(size_t)d]; e += diff*diff;
                    }
                    worst = fmax(worst,sqrt(e));
                } else if (g.state[k] != CLOSEUP_NODE_NONE) fails++;
            } else if (u > 50 && fabs(v) < 10) {
                if (g.state[k] != CLOSEUP_NODE_NONE) fails++;
            } else if (memcmp(g.pos+3*k,truth+3*k,3*sizeof(float))) fails++;
        }
        if (test < 2 || test == 4 ? st.variational_filled != 1 || st.variational_nodes != hole || worst > .1 :
                       st.variational_filled != 0 || st.variational_refused != 1) fails++;
        fprintf(stderr,"  Cartesian infill case %d: %zu holes filled, %zu refused, withheld error %.6f vox\n",
                test,st.variational_filled,st.variational_refused,worst);
        free(truth); CloseupGrid_free(&g);
    }
    fails += cs_selftest_infill_parallel();
    fails += cs_selftest_fixed_infill();
    fails += cs_selftest_provisional();
    fails += cs_selftest_rebuild();
    fprintf(stderr,"  Cartesian infill: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails != 0;
}
