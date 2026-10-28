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
    int32_t component;
    int fill_collar;
    int n, nt, cap, dimensions;
    int *rows, *cols;
    double *values, *rhs;
    double origin[5];
} CsInfillSystem;

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
        if (s->index[at] >= 0) { active++; continue; }
        if (!cs_infill_known(s->g,at,s->fill_collar) || s->g->component[at] != s->component) return;
        for (int d = 0; d < s->dimensions; d++) {
            double p = cs_infill_coordinate(s->g,s->uv,at,d);
            if (!isfinite(p)) return;
            fixed[d] += weights[k]*(p-s->origin[d]);
        }
    }
    if (!active) return;
    for (int k = 0; k < count; k++) {
        int a = s->index[nodes[k]]; if (a < 0) continue;
        for (int d = 0; d < s->dimensions; d++) s->rhs[(size_t)d*(size_t)s->n+(size_t)a] -= weights[k]*fixed[d];
        for (int j = 0; j <= k; j++) {
            int b = s->index[nodes[j]]; if (b < 0) continue;
            if (s->nt >= s->cap) { s->nt = -1; return; }
            s->rows[s->nt] = a > b ? a : b;
            s->cols[s->nt] = a > b ? b : a;
            s->values[s->nt++] = weights[k]*weights[j];
        }
    }
}

static double cs_infill_point(const CsInfillSystem *s, const double *solution, size_t at, int d)
{
    if (s->index[at] < 0) return cs_infill_coordinate(s->g,s->uv,at,d);
    double value = solution[(size_t)d*(size_t)s->n+(size_t)s->index[at]]+s->origin[d];
    /* The gates inspect the same float world points that will be exported. */
    return d < 3 ? (double)(float)value : value;
}

static int cs_infill_solve(CloseupGrid *g, const int32_t *nodes, size_t count, int32_t *index,
                            size_t x0, size_t y0, size_t x1, size_t y1, int32_t component,
                            const double normal[3], int fill_collar, int preserve_faces,
                            int provisional, double *uv)
{
    size_t box = (x1-x0+3)*(y1-y0+3);
    if (count > CS_INFILL_MAX_NODES || box > 4*CS_INFILL_MAX_NODES) return 0;
    CsInfillSystem s = {0}; SparseFactor_T factor = NULL; double *solution = NULL;
    int kept = 0;
    const char *reason = "equations";
    size_t failed_at = 0;
    double failed_lo = 0, failed_hi = 0, failed_dot = 0, failed_det = 0;
    s.g = g; s.index = index; s.component = component; s.fill_collar = fill_collar;
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
    for (size_t k = 0; k < count; k++) index[(size_t)nodes[k]] = (int32_t)k;
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
            for (int k = 0; k < 4; k++) if (index[corners[k]] >= 0) { unknown++; missing = k; }
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
                v[k] = corners[tri[t][k]]; active += index[v[k]] >= 0;
                if (index[v[k]] < 0 && (!cs_infill_known(g,v[k],fill_collar) || g->component[v[k]] != component)) valid = 0;
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
    if (!kept) fprintf(stderr,"    Cartesian infill gate: %s box=%zu,%zu..%zu,%zu at=%zu,%zu "
                      "squared-scales=%.8g,%.8g orientation=%.8g uv-det=%.8g\n",reason,x0,y0,x1,y1,
                      failed_at%g->ni,failed_at/g->ni,failed_lo,failed_hi,failed_dot,failed_det);
    for (size_t k = 0; k < count; k++) index[(size_t)nodes[k]] = -1;
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
        if (exterior || (require_fill && !has_fill)) continue;
        st->variational_holes++;
        int result = 0;
        if (!mixed && !beyond_limit && component >= 0 && count <= CS_INFILL_MAX_NODES)
            result = cs_infill_solve(g,queue,count,index,x0,y0,x1,y1,component,normal,
                                    fill_collar,preserve_faces,provisional,uv);
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
    free(seen); free(queue); free(index); return rc;
}

int CloseupSurface_fill_holes(CloseupGrid *g, double max_distance_px, CloseupFitStats *st)
{
    if (!g) return -1;
    int rc=cs_fill_holes_pass(g,max_distance_px,st,0,0,0,0,NULL,g->nrm,NULL);
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
    fails += cs_selftest_fixed_infill();
    fails += cs_selftest_provisional();
    fails += cs_selftest_rebuild();
    fprintf(stderr,"  Cartesian infill: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails != 0;
}
