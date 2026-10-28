/* Repair inferred material without treating an earlier extrapolation as an
 * observation. This is deliberately not part of complete_fixed: the original
 * artifact remains the source of a separate DATA-preservation proof. */
#include "../remesh/intersection_cleanup.h"
typedef struct CsRebuildWitness {
    int32_t ca, cb, a, b;
} CsRebuildWitness;

static int cs_rebuild_i32_compare(const void *va, const void *vb)
{
    int32_t a = *(const int32_t *)va, b = *(const int32_t *)vb;
    return a < b ? -1 : a > b;
}

static int32_t cs_rebuild_find(int32_t *parent, int32_t a)
{
    while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
    return a;
}

static int32_t cs_rebuild_label(const int32_t *labels, size_t count, int32_t label)
{
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo+(hi-lo)/2;
        if (labels[mid] < label) lo = mid+1; else hi = mid;
    }
    return lo < count && labels[lo] == label ? (int32_t)lo : -1;
}

static int cs_rebuild_witness_compare(const void *va, const void *vb)
{
    const CsRebuildWitness *a = (const CsRebuildWitness *)va, *b = (const CsRebuildWitness *)vb;
    if (a->ca != b->ca) return a->ca < b->ca ? -1 : 1;
    if (a->cb != b->cb) return a->cb < b->cb ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : a->b > b->b;
}

static void cs_rebuild_clear(CloseupGrid *g, size_t k)
{
    g->state[k] = CLOSEUP_NODE_NONE; g->chart[k] = -1; g->component[k] = -1;
    g->support[k] = 0; g->dist[k] = INFINITY;
    for (int d = 0; d < 3; d++) {
        g->pos[3*k+(size_t)d] = NAN;
        g->tu[3*k+(size_t)d] = g->tv[3*k+(size_t)d] = g->nrm[3*k+(size_t)d] = 0;
    }
}

static int cs_rebuild_copy_grid(const CloseupGrid *g, CloseupGrid *copy)
{
    size_t nn = g->ni*g->nj;
    memset(copy,0,sizeof *copy);
    copy->ni = g->ni; copy->nj = g->nj; copy->u0 = g->u0; copy->v0 = g->v0; copy->h = g->h;
    copy->pos = (float *)malloc(3*nn*sizeof(float)); copy->nrm = (float *)malloc(3*nn*sizeof(float));
    copy->tu = (float *)malloc(3*nn*sizeof(float)); copy->tv = (float *)malloc(3*nn*sizeof(float));
    copy->dist = (float *)malloc(nn*sizeof(float)); copy->support = (float *)malloc(nn*sizeof(float));
    copy->chart = (int32_t *)malloc(nn*sizeof(int32_t)); copy->component = (int32_t *)malloc(nn*sizeof(int32_t));
    copy->state = (uint8_t *)malloc(nn);
    if (!copy->pos || !copy->nrm || !copy->tu || !copy->tv || !copy->dist || !copy->support ||
        !copy->chart || !copy->component || !copy->state) { CloseupGrid_free(copy); return -1; }
    memcpy(copy->pos,g->pos,3*nn*sizeof(float)); memcpy(copy->nrm,g->nrm,3*nn*sizeof(float));
    memcpy(copy->tu,g->tu,3*nn*sizeof(float)); memcpy(copy->tv,g->tv,3*nn*sizeof(float));
    memcpy(copy->dist,g->dist,nn*sizeof(float)); memcpy(copy->support,g->support,nn*sizeof(float));
    memcpy(copy->chart,g->chart,nn*sizeof(int32_t)); memcpy(copy->component,g->component,nn*sizeof(int32_t));
    memcpy(copy->state,g->state,nn); return 0;
}

/* Four possible triangles in a regular cell, each bit names its missing
 * corner: BCD, ACD, ABD, ABC. Internal triangles have positive frame area. */
static const int cs_rebuild_triangles[4][3] = {{1,3,2},{0,3,2},{0,1,3},{0,1,2}};

static int cs_rebuild_source_faces(const CloseupGrid *g, const int32_t *faces, size_t nf,
                                    uint8_t *mask, double *tu, double *tv, float *reference,
                                    size_t *weight, int8_t *winding)
{
    for (size_t f = 0; f < nf; f++) {
        size_t v[3], x0 = g->ni, x1 = 0, y0 = g->nj, y1 = 0;
        unsigned corners = 0;
        for (int k = 0; k < 3; k++) {
            int32_t node = faces[3*f+(size_t)k];
            if (node < 0 || (size_t)node >= g->ni*g->nj || g->state[node] != CLOSEUP_NODE_DATA) return -1;
            v[k] = (size_t)node;
            size_t x = v[k]%g->ni, y = v[k]/g->ni;
            if (x < x0) x0 = x; if (x > x1) x1 = x;
            if (y < y0) y0 = y; if (y > y1) y1 = y;
        }
        if (x1 != x0+1 || y1 != y0+1 || g->component[v[0]] != g->component[v[1]] ||
            g->component[v[0]] != g->component[v[2]]) return -1;
        for (int k = 0; k < 3; k++) corners |= 1u << (unsigned)((v[k]/g->ni-y0)*2+v[k]%g->ni-x0);
        unsigned missing = 0;
        for (; missing < 4 && (corners & (1u << missing)); missing++) { }
        if (missing == 4 || corners != (15u^(1u << missing))) return -1;
        size_t cell = y0*g->ni+x0;
        if (mask[cell] & (1u << missing)) return -1;
        /* Only the two complementary halves of one diagonal may coexist. */
        unsigned diagonal = missing == 0 || missing == 3 ? 9u : 6u;
        if (mask[cell] & ~diagonal) return -1;
        mask[cell] |= (uint8_t)(1u << missing);
        double au = g->h*((double)(v[1]%g->ni)-(double)(v[0]%g->ni));
        double av = g->h*((double)(v[1]/g->ni)-(double)(v[0]/g->ni));
        double bu = g->h*((double)(v[2]%g->ni)-(double)(v[0]%g->ni));
        double bv = g->h*((double)(v[2]/g->ni)-(double)(v[0]/g->ni));
        double det = au*bv-av*bu, ut[3], vt[3];
        int8_t sign = det > 0 ? 1 : -1;
        int32_t component = g->component[v[0]];
        if (winding[component] && winding[component] != sign) return -1;
        winding[component] = sign;
        for (int d = 0; d < 3; d++) {
            double e = (double)g->pos[3*v[1]+(size_t)d]-g->pos[3*v[0]+(size_t)d];
            double a = (double)g->pos[3*v[2]+(size_t)d]-g->pos[3*v[0]+(size_t)d];
            ut[d] = (e*bv-a*av)/det; vt[d] = (a*au-e*bu)/det;
        }
        for (int k = 0; k < 3; k++) {
            for (int d = 0; d < 3; d++) { tu[3*v[k]+(size_t)d] += ut[d]; tv[3*v[k]+(size_t)d] += vt[d]; }
            weight[v[k]]++;
        }
    }
    for (size_t k = 0; k < g->ni*g->nj; k++) if (weight[k]) {
        double normal[3], length = 0;
        for (int d = 0; d < 3; d++) { tu[3*k+(size_t)d] /= (double)weight[k]; tv[3*k+(size_t)d] /= (double)weight[k]; }
        for (int d = 0; d < 3; d++) {
            normal[d] = tu[3*k+(size_t)((d+1)%3)]*tv[3*k+(size_t)((d+2)%3)]-
                        tu[3*k+(size_t)((d+2)%3)]*tv[3*k+(size_t)((d+1)%3)];
            length += normal[d]*normal[d];
        }
        if (!(length > 0) || !isfinite(length)) return -1;
        length = sqrt(length);
        for (int d = 0; d < 3; d++) reference[3*k+(size_t)d] = (float)(normal[d]/length);
    }
    return 0;
}

static int cs_rebuild_witness_pass(const CloseupGrid *g, size_t a, size_t b,
                                   const double *tu, const double *tv, const float *normal,
                                   const size_t *weight, const int8_t *winding)
{
    if (!weight[a] || !weight[b] || winding[g->component[a]] != winding[g->component[b]]) return 0;
    double du = g->h*((double)(b%g->ni)-(double)(a%g->ni));
    double dv = g->h*((double)(b/g->ni)-(double)(a/g->ni));
    double frame = hypot(du,dv), length2 = 0, residual2 = 0, oa = 0, ob = 0, dot = 0;
    for (int d = 0; d < 3; d++) {
        double delta = (double)g->pos[3*b+(size_t)d]-g->pos[3*a+(size_t)d];
        double prediction = .5*((tu[3*a+(size_t)d]+tu[3*b+(size_t)d])*du+
                                 (tv[3*a+(size_t)d]+tv[3*b+(size_t)d])*dv);
        length2 += delta*delta; residual2 += (delta-prediction)*(delta-prediction);
        oa += normal[3*a+(size_t)d]*delta; ob += normal[3*b+(size_t)d]*delta;
        dot += (double)normal[3*a+(size_t)d]*normal[3*b+(size_t)d];
    }
    return sqrt(length2) <= CLOSEUP_BRIDGE_STRETCH*frame+2 &&
           fabs(oa) <= CLOSEUP_BRIDGE_NORMAL_VOX && fabs(ob) <= CLOSEUP_BRIDGE_NORMAL_VOX &&
           sqrt(residual2) <= CLOSEUP_BRIDGE_NORMAL_VOX && dot >= CLOSEUP_DEDUP_COS;
}

/* Packed-frame proximity only proposes a witness; chart identity and every
 * measured physical continuation must agree before a component union. A
 * single point pair is insufficient evidence for an oriented boundary. */
static int cs_rebuild_groups(const CloseupGrid *g, const uint8_t *domain, const int32_t *nearest,
                              const double *tu, const double *tv, const float *normal,
                              const size_t *weight, const int8_t *winding, int32_t *parent,
                              CloseupRebuildStats *st)
{
    CsRebuildWitness *witness = NULL;
    size_t count = 0, cap = 0;
    int rc = -1;
    const int dx[4] = {1,0,1,-1}, dy[4] = {0,1,1,1};
    double max_gap2 = 4*CLOSEUP_BRIDGE_PX*CLOSEUP_BRIDGE_PX;
    for (size_t k = 0; k < g->ni*g->nj; k++) if (domain[k] && nearest[k] >= 0) {
        size_t x = k%g->ni, y = k/g->ni;
        for (int direction = 0; direction < 4; direction++) {
            if ((dx[direction] < 0 && !x) || (dx[direction] > 0 && x+1 >= g->ni) ||
                (dy[direction] > 0 && y+1 >= g->nj)) continue;
            size_t q = (size_t)((int64_t)k+(int64_t)dy[direction]*(int64_t)g->ni+dx[direction]);
            if (!domain[q] || nearest[q] < 0) continue;
            int32_t a = nearest[k], b = nearest[q], ca = g->component[a], cb = g->component[b];
            if (ca == cb || g->chart[a] < 0 || g->chart[a] != g->chart[b]) continue;
            double du = g->h*((double)((size_t)a%g->ni)-(double)((size_t)b%g->ni));
            double dv = g->h*((double)((size_t)a/g->ni)-(double)((size_t)b/g->ni));
            if (du*du+dv*dv > max_gap2) continue;
            if (ca > cb) { int32_t swap = ca; ca = cb; cb = swap; swap = a; a = b; b = swap; }
            if (count == cap) {
                size_t next = cap ? cap*2 : 256;
                if (next < cap || next > SIZE_MAX/sizeof *witness) goto done;
                CsRebuildWitness *grown = (CsRebuildWitness *)realloc(witness,next*sizeof *witness);
                if (!grown) goto done;
                witness = grown; cap = next;
            }
            witness[count].ca = ca; witness[count].cb = cb; witness[count].a = a; witness[count].b = b; count++;
        }
    }
    if (count) qsort(witness,count,sizeof *witness,cs_rebuild_witness_compare);
    for (size_t begin = 0; begin < count;) {
        size_t end = begin, unique = 0, passed = 0;
        int distinct_a = 0, distinct_b = 0;
        while (end < count && witness[end].ca == witness[begin].ca && witness[end].cb == witness[begin].cb) end++;
        for (size_t k = begin; k < end; k++) {
            if (k > begin && !cs_rebuild_witness_compare(witness+k,witness+k-1)) continue;
            unique++;
            distinct_a |= witness[k].a != witness[begin].a; distinct_b |= witness[k].b != witness[begin].b;
            passed += (size_t)cs_rebuild_witness_pass(g,(size_t)witness[k].a,(size_t)witness[k].b,tu,tv,normal,weight,winding);
        }
        int qualified = unique >= 2 && distinct_a && distinct_b && passed == unique;
        if (qualified) {
            int32_t a = cs_rebuild_find(parent,witness[begin].ca), b = cs_rebuild_find(parent,witness[begin].cb);
            if (a != b) { parent[b] = a; st->certified_joins++; }
        } else st->refused_joins++;
        fprintf(stderr,"  DATA rebuild component witness %d,%d: %zu/%zu physical pairs, %s\n",
                witness[begin].ca,witness[begin].cb,passed,unique,qualified ? "certified" : "separate");
        begin = end;
    }
    rc = 0;
done:
    free(witness); return rc;
}

/* Return the exact non-overlapping cell triangles allowed by its observed
 * diagonal and present corners. The caller omits observed triangles. */
static int cs_rebuild_cell_triangles(const CloseupGrid *g, size_t cell, uint8_t observed, int out[2])
{
    size_t corners[4] = {cell,cell+1,cell+g->ni,cell+g->ni+1};
    int valid = 0, missing = -1;
    for (int k = 0; k < 4; k++) { valid += g->state[corners[k]] != 0; if (!g->state[corners[k]]) missing = k; }
    if (valid < 3) return 0;
    if (valid == 3) { out[0] = missing; return 1; }
    if (observed & 9u) { out[0] = 3; out[1] = 0; }
    else { out[0] = 2; out[1] = 1; }
    return 2;
}

/* Natural boundary: second-derivative stencils crossing the bounded material
 * domain or a separated component are omitted. A DATA-tangent screening term
 * fixes the unconstrained affine modes of sparse/one-sided domains. Its smooth
 * inverse-distance-fourth weight has the units of the Hessian energy, and
 * supplies no old generated geometry. No mesh-density threshold is used. */
static int cs_rebuild_solve(CloseupGrid *g, const int32_t *nodes, size_t count, int32_t *index,
                             size_t x0, size_t y0, size_t x1, size_t y1, int32_t component,
                             const int32_t *nearest, const double *tu, const double *tv,
                             CloseupRebuildStats *st)
{
    size_t left = x0 ? x0-1 : 0, top = y0 ? y0-1 : 0;
    size_t right = x1+1 < g->ni ? x1+1 : x1, bottom = y1+1 < g->nj ? y1+1 : y1;
    size_t box = (right-left+1)*(bottom-top+1);
    CsInfillSystem s = {0}; SparseFactor_T factor = NULL; double *solution = NULL;
    int rc = -1;
    if (count > INT_MAX || box > ((size_t)INT_MAX-count)/22) return -1;
    s.g = g; s.index = index; s.component = component; s.dimensions = 3;
    s.n = (int)count; s.cap = (int)(22*box+count);
    s.rows = (int *)malloc((size_t)s.cap*sizeof(int)); s.cols = (int *)malloc((size_t)s.cap*sizeof(int));
    s.values = (double *)malloc((size_t)s.cap*sizeof(double)); s.rhs = (double *)calloc(3*count,sizeof(double));
    solution = (double *)malloc(3*count*sizeof(double));
    if (!s.rows || !s.cols || !s.values || !s.rhs || !solution) goto done;
    if (!count || nearest[nodes[0]] < 0) { rc = 0; goto done; }
    for (int d = 0; d < 3; d++) s.origin[d] = g->pos[3*(size_t)nearest[nodes[0]]+(size_t)d];
    for (size_t k = 0; k < count; k++) index[nodes[k]] = (int32_t)k;
    for (size_t j = top; j <= bottom; j++) for (size_t i = left; i <= right; i++) {
        size_t at = j*g->ni+i;
        const double second[3] = {1,-2,1};
        const double mixed[4] = {.3535533905932738,-.3535533905932738,-.3535533905932738,.3535533905932738};
        if (i && i+1 < g->ni) { size_t row[3] = {at-1,at,at+1}; cs_infill_row(&s,row,second,3); }
        if (j && j+1 < g->nj) { size_t row[3] = {at-g->ni,at,at+g->ni}; cs_infill_row(&s,row,second,3); }
        if (i && j && i+1 < g->ni && j+1 < g->nj) {
            size_t row[4] = {at-g->ni-1,at-g->ni+1,at+g->ni-1,at+g->ni+1};
            cs_infill_row(&s,row,mixed,4);
        }
        if (s.nt < 0) goto done;
    }
    for (size_t k = 0; k < count; k++) {
        size_t at = (size_t)nodes[k]; int32_t source = nearest[at];
        if (source < 0 || g->state[source] != CLOSEUP_NODE_DATA || g->component[source] != component) { rc = 0; goto done; }
        double du = g->h*((double)(at%g->ni)-(double)((size_t)source%g->ni));
        double dv = g->h*((double)(at/g->ni)-(double)((size_t)source/g->ni));
        double ratio = g->h*g->h/(g->h*g->h+du*du+dv*dv), lambda = ratio*ratio;
        if (!(lambda > 0) || !isfinite(lambda)) { rc = 0; goto done; }
        s.rows[s.nt] = s.cols[s.nt] = (int)k; s.values[s.nt++] = lambda;
        for (int d = 0; d < 3; d++) {
            double target = (double)g->pos[3*(size_t)source+(size_t)d]+tu[3*(size_t)source+(size_t)d]*du+
                            tv[3*(size_t)source+(size_t)d]*dv-s.origin[d];
            if (!isfinite(target)) { rc = 0; goto done; }
            s.rhs[(size_t)d*count+k] += lambda*target;
        }
    }
    fprintf(stderr,"  DATA rebuild screened solve: %zu unknowns, %d triplets, group %d, box=%zu,%zu..%zu,%zu\n",
            count,s.nt,component,x0,y0,x1,y1);
    /* Allocation/factor failures are retryable runtime errors, not geometric
     * refusals or proof that a large region cannot be repaired. */
    if (Sparse_factor_spd(s.n,s.nt,s.rows,s.cols,s.values,&factor) ||
        Sparse_factor_solve_multi(factor,s.rhs,solution,3)) goto done;
    for (size_t k = 0; k < count; k++) for (int d = 0; d < 3; d++) {
        double value = cs_infill_point(&s,solution,(size_t)nodes[k],d);
        if (!isfinite(value)) { rc = 0; goto done; }
    }
    for (size_t k = 0; k < count; k++) {
        size_t at = (size_t)nodes[k];
        size_t source = (size_t)nearest[at];
        double du = g->h*((double)(at%g->ni)-(double)(source%g->ni));
        double dv = g->h*((double)(at/g->ni)-(double)(source/g->ni));
        double ratio = g->h*g->h/(g->h*g->h+du*du+dv*dv);
        for (int d = 0; d < 3; d++) {
            g->pos[3*at+(size_t)d] = (float)cs_infill_point(&s,solution,at,d);
            double target = (double)g->pos[3*source+(size_t)d]+tu[3*source+(size_t)d]*du+tv[3*source+(size_t)d]*dv;
            double error = (double)g->pos[3*at+(size_t)d]-target;
            st->tangent_screen_energy += ratio*ratio*error*error;
        }
        g->state[at] = CLOSEUP_NODE_FILL; g->chart[at] = -1; g->support[at] = 0;
    }
    st->screened_nodes += count;
    rc = 1;
done:
    for (size_t k = 0; k < count; k++) index[nodes[k]] = -1;
    Sparse_factor_free(&factor); free(s.rows); free(s.cols); free(s.values); free(s.rhs); free(solution);
    return rc;
}

/* A local DATA orientation reference is essential on rolled sheets: the mean
 * normal over a large patch can vanish or oppose valid triangles. */
static int cs_rebuild_triangle_oriented(const CloseupGrid *g, const size_t v[3],
                                         const int32_t *nearest, const float *reference, double diagnostic[3])
{
    double e[3], f[3], normal[3], length2 = 0, dot = 0, ref[3] = {0};
    double source_normal[3][3] = {{0}}, coherence = 1, ref_length2 = 0;
    if (diagnostic) diagnostic[0] = diagnostic[1] = diagnostic[2] = 0;
    for (int k = 0; k < 3; k++) {
        int32_t data = nearest[v[k]];
        if (data < 0) return 0;
        double reference_length2 = 0;
        for (int d = 0; d < 3; d++) {
            double n = reference[3*(size_t)data+(size_t)d];
            reference_length2 += n*n; ref[d] += n; source_normal[k][d] = n;
        }
        if (!(reference_length2 > 0)) return 0;
    }
    for (int d = 0; d < 3; d++) {
        e[d] = (double)g->pos[3*v[1]+(size_t)d]-g->pos[3*v[0]+(size_t)d];
        f[d] = (double)g->pos[3*v[2]+(size_t)d]-g->pos[3*v[0]+(size_t)d];
    }
    for (int d = 0; d < 3; d++) {
        normal[d] = e[(d+1)%3]*f[(d+2)%3]-e[(d+2)%3]*f[(d+1)%3];
        length2 += normal[d]*normal[d]; dot += normal[d]*ref[d]; ref_length2 += ref[d]*ref[d];
    }
    for (int a = 0; a < 3; a++) for (int b = a+1; b < 3; b++) {
        double similarity = 0;
        for (int d = 0; d < 3; d++) similarity += source_normal[a][d]*source_normal[b][d];
        if (similarity < coherence) coherence = similarity;
    }
    if (diagnostic) {
        diagnostic[0] = length2;
        diagnostic[1] = length2 > 0 && ref_length2 > 0 ? dot/sqrt(length2*ref_length2) : 0;
        diagnostic[2] = coherence;
    }
    return isfinite(length2) && length2 > 0 && isfinite(dot) && dot > 0;
}

/* The observed diagonal is immutable. A wholly inferred cell may use either
 * triangulation: choose a feasible one, then the smaller worst triangle shape
 * ratio. This changes no source point and does not waive any orientation gate. */
static int cs_rebuild_choose_diagonal(const CloseupGrid *g, size_t cell, uint8_t observed,
                                       const int32_t *nearest, const float *reference,
                                       int triangles[2], CloseupRebuildStats *st)
{
    int count = cs_rebuild_cell_triangles(g,cell,observed,triangles);
    if (observed || count != 2) return count;
    size_t corners[4] = {cell,cell+1,cell+g->ni,cell+g->ni+1};
    int32_t component = g->component[cell];
    for (int k = 1; k < 4; k++) if (g->component[corners[k]] != component) return count;
    const int options[2][2] = {{2,1},{3,0}};
    int feasible[2] = {1,1}; double worst[2] = {0,0};
    int32_t local_faces[2][2][3] = {{{0}}};
    (void)nearest; (void)reference;
    for (int option = 0; option < 2; option++) for (int t = 0; t < 2; t++) {
        size_t v[3]; double area2 = 0, edges2 = 0, a[3], b[3];
        for (int k = 0; k < 3; k++) {
            v[k] = corners[cs_rebuild_triangles[options[option][t]][k]];
            local_faces[option][t][k] = (int32_t)v[k];
        }
        for (int d = 0; d < 3; d++) {
            a[d] = (double)g->pos[3*v[1]+(size_t)d]-g->pos[3*v[0]+(size_t)d];
            b[d] = (double)g->pos[3*v[2]+(size_t)d]-g->pos[3*v[0]+(size_t)d];
        }
        for (int d = 0; d < 3; d++) {
            double value = a[(d+1)%3]*b[(d+2)%3]-a[(d+2)%3]*b[(d+1)%3];
            area2 += value*value;
        }
        if (!(area2 > 0) || !isfinite(area2)) { feasible[option] = 0; continue; }
        for (int k = 0; k < 3; k++) for (int d = 0; d < 3; d++) {
            double delta = (double)g->pos[3*v[(k+1)%3]+(size_t)d]-g->pos[3*v[k]+(size_t)d];
            edges2 += delta*delta;
        }
        double shape = edges2/sqrt(area2);
        if (shape > worst[option]) worst[option] = shape;
    }
    IntersectionCleanupParams params; IntersectionCleanup_default_params(&params);
    for (int option = 0; option < 2; option++) {
        if (IntersectionCleanup_pair_test(g->pos,g->ni*g->nj,local_faces[option][0],local_faces[option][1],&params,NULL))
            feasible[option] = 0;
    }
    if (feasible[1] && (!feasible[0] || worst[1] < worst[0])) {
        triangles[0] = 3; triangles[1] = 0; st->alternate_diagonals++;
    }
    return count;
}

/* All inferred triangles sharing an inferred edge form one atomic patch.
 * This also joins an all-DATA cap to its neighboring generated geometry. */
static int cs_rebuild_patch_components(const CloseupGrid *g, CloseupRebuildResult *out,
                                        uint8_t *refused, size_t *npatches)
{
    size_t nn = g->ni*g->nj, n = *npatches, groups = 0;
    int32_t *edge = (int32_t *)malloc(4*nn*sizeof(int32_t));
    int32_t *parent = (int32_t *)malloc((n ? n : 1)*sizeof(int32_t));
    int32_t *map = (int32_t *)malloc((n ? n : 1)*sizeof(int32_t));
    uint8_t *state = (uint8_t *)calloc(n ? n : 1,1);
    int rc = -1;
    if (!edge || !parent || !map || !state) goto done;
    for (size_t k = 0; k < 4*nn; k++) edge[k] = -1;
    for (size_t p = 0; p < n; p++) { parent[p] = (int32_t)p; map[p] = -1; }
    for (size_t f = 0; f < out->faces; f++) for (int k = 0; k < 3; k++) {
        size_t a = (size_t)out->face_nodes[3*f+(size_t)k], b = (size_t)out->face_nodes[3*f+(size_t)((k+1)%3)];
        if (a > b) { size_t swap = a; a = b; b = swap; }
        size_t delta = b-a, kind = delta == 1 ? 0 : delta == g->ni ? 1 : delta == g->ni+1 ? 2 : 3;
        size_t key = 4*a+kind;
        int32_t p = out->face_patch[f];
        if (edge[key] < 0) edge[key] = p;
        else {
            int32_t pa = cs_rebuild_find(parent,p), pb = cs_rebuild_find(parent,edge[key]);
            if (pa != pb) parent[pb] = pa;
        }
    }
    for (size_t p = 0; p < n; p++) {
        int32_t root = cs_rebuild_find(parent,(int32_t)p);
        state[root] |= refused[p];
        if (map[root] < 0) map[root] = (int32_t)groups++;
    }
    memset(refused,0,n);
    for (size_t p = 0; p < n; p++) if (parent[p] == (int32_t)p) refused[map[p]] = state[p];
    for (size_t k = 0; k < nn; k++) if (out->node_patch[k] >= 0)
        out->node_patch[k] = map[cs_rebuild_find(parent,out->node_patch[k])];
    for (size_t f = 0; f < out->faces; f++) out->face_patch[f] = map[cs_rebuild_find(parent,out->face_patch[f])];
    *npatches = groups; rc = 0;
done:
    free(edge); free(parent); free(map); free(state); return rc;
}

#include "closeup_rebuild_motion.inc"
#include "closeup_rebuild_transport.inc"
#include "closeup_rebuild_domain.inc"

void CloseupRebuildResult_free(CloseupRebuildResult *result)
{
    if (!result) return;
    free(result->node_patch); free(result->component_group); free(result->face_nodes); free(result->face_patch);
    memset(result,0,sizeof *result);
}

int CloseupSurface_rebuild_generated(CloseupGrid *g, const uint8_t *domain,
                                     const int32_t *observed_faces, size_t observed_nf,
                                     double max_distance_px, CloseupRebuildResult *result,
                                     CloseupRebuildStats *st)
{
    CloseupGrid work = {0}, audit = {0}; CloseupRebuildResult out = {0};
    uint8_t *seed = NULL, *mask = NULL, *seen = NULL, *refused = NULL;
    float *d2 = NULL, *reference = NULL, *target = NULL;
    double *tu = NULL, *tv = NULL;
    size_t *weight = NULL;
    int32_t *nearest = NULL, *labels = NULL, *parent = NULL, *queue = NULL, *index = NULL, *remap = NULL;
    int8_t *winding = NULL;
    size_t nn = 0, ncomponents = 0, nlabels = 0, npatches = 0;
    int rc = -1;
    if (!g || !domain || !result || !st || (observed_nf && !observed_faces) || !g->ni || !g->nj ||
        g->ni > INT32_MAX/g->nj || !(g->h > 0) || !isfinite(g->h) || !(max_distance_px > 0) ||
        !isfinite(max_distance_px)) return -1;
    memset(result,0,sizeof *result); memset(st,0,sizeof *st); st->metric_qualification_pending = 1;
    nn = g->ni*g->nj;
    if (nn > SIZE_MAX/(22*sizeof(double)) || nn > (size_t)INT_MAX/3) return -1;
    seed = (uint8_t *)calloc(nn,1); mask = (uint8_t *)calloc(nn,1); seen = (uint8_t *)calloc(nn,1);
    d2 = (float *)malloc(nn*sizeof(float)); reference = (float *)calloc(3*nn,sizeof(float));
    tu = (double *)calloc(3*nn,sizeof(double)); tv = (double *)calloc(3*nn,sizeof(double));
    weight = (size_t *)calloc(nn,sizeof(size_t)); nearest = (int32_t *)malloc(nn*sizeof(int32_t));
    labels = (int32_t *)malloc(nn*sizeof(int32_t)); parent = (int32_t *)malloc(nn*sizeof(int32_t));
    queue = (int32_t *)malloc(nn*sizeof(int32_t)); index = (int32_t *)malloc(nn*sizeof(int32_t));
    winding = (int8_t *)calloc(nn,1);
    out.node_patch = (int32_t *)malloc(nn*sizeof(int32_t)); out.component_group = (int32_t *)malloc(nn*sizeof(int32_t));
    if (!seed || !mask || !seen || !d2 || !reference || !tu || !tv || !weight || !nearest || !labels ||
        !parent || !queue || !index || !winding || !out.node_patch || !out.component_group || cs_rebuild_copy_grid(g,&work)) goto done;
    for (size_t k = 0; k < nn; k++) {
        index[k] = out.node_patch[k] = out.component_group[k] = -1;
        st->domain_nodes += domain[k] != 0;
        st->enclosed_gap_nodes += domain[k] && !g->state[k];
        if (g->state[k] && (!domain[k] || (g->state[k] != CLOSEUP_NODE_DATA && g->state[k] != CLOSEUP_NODE_FILL))) goto done;
        if (g->state[k] == CLOSEUP_NODE_DATA) {
            if (g->component[k] < 0) goto done;
            for (int d = 0; d < 3; d++) if (!isfinite(g->pos[3*k+(size_t)d]) || !isfinite(g->nrm[3*k+(size_t)d])) goto done;
            seed[k] = 1; labels[nlabels++] = g->component[k]; st->data_nodes++;
        } else {
            st->old_fill_nodes += g->state[k] == CLOSEUP_NODE_FILL; cs_rebuild_clear(&work,k);
        }
    }
    if (nlabels) qsort(labels,nlabels,sizeof *labels,cs_rebuild_i32_compare);
    for (size_t k = 0; k < nlabels; k++) if (!k || labels[k] != labels[k-1]) labels[ncomponents++] = labels[k];
    st->source_components = ncomponents;
    for (size_t c = 0; c < ncomponents; c++) parent[c] = (int32_t)c;
    for (size_t k = 0; k < nn; k++) if (seed[k]) work.component[k] = cs_rebuild_label(labels,ncomponents,g->component[k]);
    if (cs_edt2d(seed,g->ni,g->nj,d2,nearest) ||
        cs_rebuild_source_faces(&work,observed_faces,observed_nf,mask,tu,tv,reference,weight,winding)) goto done;
    if (cs_rebuild_groups(&work,domain,nearest,tu,tv,reference,weight,winding,parent,st)) goto done;
    for (size_t c = 0; c < ncomponents; c++) st->supported_groups += cs_rebuild_find(parent,(int32_t)c) == (int32_t)c;
    for (size_t k = 0; k < nn; k++) if (domain[k] && nearest[k] >= 0) {
        int32_t c = cs_rebuild_find(parent,work.component[nearest[k]]);
        out.component_group[k] = c;
    }
    for (size_t k = 0; k < nn; k++) if (domain[k]) {
        work.component[k] = out.component_group[k];
        if (!seed[k]) work.dist[k] = nearest[k] >= 0 ? (float)(g->h*sqrt((double)d2[k])) : INFINITY;
    }
    for (size_t start = 0; start < nn; start++) {
        if (!domain[start] || seed[start] || seen[start]) continue;
        size_t count = 1, x0 = g->ni, x1 = 0, y0 = g->nj, y1 = 0;
        int unsupported = 0, beyond = 0;
        int32_t component = work.component[start];
        seen[start] = 1; queue[0] = (int32_t)start;
        for (size_t head = 0; head < count; head++) {
            size_t k = (size_t)queue[head], x = k%g->ni, y = k/g->ni;
            if (x < x0) x0 = x; if (x > x1) x1 = x;
            if (y < y0) y0 = y; if (y > y1) y1 = y;
            beyond |= !(work.dist[k] <= max_distance_px);
            unsupported |= component < 0 || nearest[k] < 0 || !weight[nearest[k]];
            for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
                if ((!dx && !dy) || (dx < 0 && !x) || (dy < 0 && !y) ||
                    (dx > 0 && x+1 >= g->ni) || (dy > 0 && y+1 >= g->nj)) continue;
                size_t q = (size_t)((int64_t)k+(int64_t)dy*(int64_t)g->ni+dx);
                if (!domain[q] || seed[q] || seen[q] || work.component[q] != component) continue;
                seen[q] = 1; queue[count++] = (int32_t)q;
            }
        }
        st->patches++;
        int solved = 0;
        if (!unsupported && !beyond) solved = cs_rebuild_solve(&work,queue,count,index,x0,y0,x1,y1,component,nearest,tu,tv,st);
        if (solved < 0) goto done;
        if (!solved) {
            st->refused_patches++; st->reach_refused += beyond != 0; st->unsupported_nodes += unsupported ? count : 0;
            for (size_t k = 0; k < count; k++) cs_rebuild_clear(&work,(size_t)queue[k]);
            fprintf(stderr,"  DATA rebuild refused %zu nodes: unsupported=%d beyond-reach=%d\n",count,unsupported,beyond);
            continue;
        }
        for (size_t k = 0; k < count; k++) out.node_patch[queue[k]] = (int32_t)npatches;
        npatches++;
    }
    /* Stage the DATA-only variational result. The historical distant-normal
     * motion is retained above for regression evidence, but must not distort
     * the sheet to agree with a donor normal far along a curved surface. */
    target = (float *)malloc(3*nn*sizeof(float));
    if (!target) goto done;
    memcpy(target,work.pos,3*nn*sizeof(float));
    /* The maximum face count follows grid representation, not a policy cap. */
    out.face_nodes = (int32_t *)malloc(6*nn*sizeof(int32_t));
    out.face_patch = (int32_t *)malloc(2*nn*sizeof(int32_t));
    refused = (uint8_t *)calloc(3*nn,1); remap = (int32_t *)malloc(3*nn*sizeof(int32_t));
    if (!out.face_nodes || !out.face_patch || !refused || !remap) goto done;
    for (size_t j = 0; j+1 < g->nj; j++) for (size_t i = 0; i+1 < g->ni; i++) {
        size_t cell = j*g->ni+i, corners[4] = {cell,cell+1,cell+g->ni,cell+g->ni+1};
        int triangles[2] = {0}, nt = cs_rebuild_choose_diagonal(&work,cell,mask[cell],nearest,reference,triangles,st), separated = 0;
        for (int t = 0; t < nt; t++) {
            int triangle = triangles[t];
            if (mask[cell] & (1u << (unsigned)triangle)) continue;
            size_t v[3]; int32_t patch = -1, component = -1; int eligible = 1, patch_conflict = 0;
            for (int k = 0; k < 3; k++) {
                v[k] = corners[cs_rebuild_triangles[triangle][k]];
                int32_t c = work.component[v[k]], p = out.node_patch[v[k]];
                if (!domain[v[k]] || !work.state[v[k]] || c < 0 || (component >= 0 && c != component)) eligible = 0;
                if (patch >= 0 && p >= 0 && p != patch) patch_conflict = 1;
                if (p >= 0) patch = p;
                if (component < 0) component = c;
            }
            if (!eligible) { separated = 1; continue; }
            if (patch_conflict) goto done;
            if (patch < 0) { patch = (int32_t)npatches++; st->patches++; }
            for (int k = 0; k < 3; k++) out.face_nodes[3*out.faces+(size_t)k] = (int32_t)v[k];
            if (winding[component] < 0) {
                int32_t swap = out.face_nodes[3*out.faces+1];
                out.face_nodes[3*out.faces+1] = out.face_nodes[3*out.faces+2]; out.face_nodes[3*out.faces+2] = swap;
            }
            out.face_patch[out.faces++] = patch;
        }
        st->separated_cells += separated != 0;
    }
    if (cs_rebuild_patch_components(&work,&out,refused,&npatches)) goto done;
    st->patches = st->refused_patches+npatches;
    CsRebuildTransport transport = {0};
    st->target_start_patches = npatches;
    if (cs_rebuild_transport_certificate(&work,&out,observed_faces,observed_nf,refused,&transport)) goto done;
    st->transport_components = transport.components; st->transport_anchored_components = transport.anchored_components;
    st->transport_unanchored_components = transport.unanchored_components; st->transport_unanchored_faces = transport.unanchored_faces;
    st->transport_degenerate_faces = transport.degenerate_faces; st->transport_bad_edges = transport.bad_edges;
    st->transport_winding_edges = transport.winding_edges; st->transport_nonmanifold_edges = transport.nonmanifold_edges;
    st->world_exact = transport.exact; st->world_proximity = transport.proximity;
    fprintf(stderr,"  DATA rebuild local transport: %zu components (%zu observed-anchored, %zu unproven); "
            "%zu unproven inferred faces, %zu degenerate faces, %zu nonpositive normal edges, %zu winding conflicts, %zu nonmanifold edges\n",
            transport.components,transport.anchored_components,transport.unanchored_components,transport.unanchored_faces,
            transport.degenerate_faces,transport.bad_edges,transport.winding_edges,transport.nonmanifold_edges);
    fprintf(stderr,"  DATA rebuild complete world contact: exact %zu old/old + %zu old/new + %zu new/new; "
            "existing proximity %zu old/old + %zu old/new + %zu new/new (original observed faces retained)\n",
            transport.exact.observed_observed,transport.exact.observed_inferred,transport.exact.inferred_inferred,
            transport.proximity.observed_observed,transport.proximity.observed_inferred,transport.proximity.inferred_inferred);
    /* New normals come from actual proposed triangles, signed like nearby
     * saved source normals. Original shading normals are never recomputed. */
    memset(tu,0,3*nn*sizeof(double)); memset(tv,0,3*nn*sizeof(double)); memset(weight,0,nn*sizeof(size_t));
    for (size_t f = 0; f < out.faces; f++) {
        size_t v[3] = {(size_t)out.face_nodes[3*f],(size_t)out.face_nodes[3*f+1],(size_t)out.face_nodes[3*f+2]};
        double au = g->h*((double)(v[1]%g->ni)-(double)(v[0]%g->ni));
        double av = g->h*((double)(v[1]/g->ni)-(double)(v[0]/g->ni));
        double bu = g->h*((double)(v[2]%g->ni)-(double)(v[0]%g->ni));
        double bv = g->h*((double)(v[2]/g->ni)-(double)(v[0]/g->ni));
        double det = au*bv-av*bu;
        for (int k = 0; k < 3; k++) if (!seed[v[k]]) {
            weight[v[k]]++;
            for (int d = 0; d < 3; d++) {
                double e = (double)work.pos[3*v[1]+(size_t)d]-work.pos[3*v[0]+(size_t)d];
                double a = (double)work.pos[3*v[2]+(size_t)d]-work.pos[3*v[0]+(size_t)d];
                tu[3*v[k]+(size_t)d] += (e*bv-a*av)/det; tv[3*v[k]+(size_t)d] += (a*au-e*bu)/det;
            }
        }
    }
    for (size_t k = 0; k < nn; k++) if (out.node_patch[k] >= 0) {
        if (!weight[k]) { refused[out.node_patch[k]] = 1; continue; }
        for (int d = 0; d < 3; d++) {
            work.tu[3*k+(size_t)d] = (float)(tu[3*k+(size_t)d]/(double)weight[k]);
            work.tv[3*k+(size_t)d] = (float)(tv[3*k+(size_t)d]/(double)weight[k]);
        }
        cs_set_normal(&work,k);
        size_t data = (size_t)nearest[k]; double sign = 0;
        for (int d = 0; d < 3; d++) sign += (double)g->nrm[3*data+(size_t)d]*reference[3*data+(size_t)d];
        if (sign < 0) for (int d = 0; d < 3; d++) work.nrm[3*k+(size_t)d] = -work.nrm[3*k+(size_t)d];
    }
    for (size_t p = 0; p < npatches; p++) if (refused[p] == 1) st->orientation_refused++;
    /* Contact sees only retained DATA and eligible proposals. Rejected
     * geometry cannot block a different patch that survived orientation. */
    if (cs_rebuild_copy_grid(&work,&audit)) goto done;
    for (size_t k = 0; k < nn; k++) if (out.node_patch[k] >= 0 && refused[out.node_patch[k]]) cs_rebuild_clear(&audit,k);
    CloseupDedupStats contact = {0};
    if (CloseupSurface_dedup_fill(&audit,CLOSEUP_DEDUP_VOX,CLOSEUP_DEDUP_FRAME_PX,CLOSEUP_DEDUP_COS,&contact)) goto done;
    for (size_t k = 0; k < nn; k++) if (out.node_patch[k] >= 0 && !refused[out.node_patch[k]] && !audit.state[k]) refused[out.node_patch[k]] = 2;
    for (size_t p = 0; p < npatches; p++) {
        if (refused[p]) {
            remap[p] = -1; st->refused_patches++; st->contact_refused += refused[p] == 2 || refused[p] == 3;
        } else remap[p] = (int32_t)out.patches++;
    }
    for (size_t k = 0; k < nn; k++) if (out.node_patch[k] >= 0) {
        out.node_patch[k] = remap[out.node_patch[k]];
        if (out.node_patch[k] < 0) cs_rebuild_clear(&work,k); else out.nodes++;
    }
    size_t kept_faces = 0;
    for (size_t f = 0; f < out.faces; f++) {
        int32_t patch = remap[out.face_patch[f]];
        if (patch < 0) continue;
        memcpy(out.face_nodes+3*kept_faces,out.face_nodes+3*f,3*sizeof(int32_t));
        out.face_patch[kept_faces++] = patch;
    }
    out.faces = kept_faces;
    /* Only now commit the disposable generated portion. DATA, including
     * chart/component/support and saved tangents, remains byte-exact. */
    for (size_t k = 0; k < nn; k++) if (!seed[k]) {
        memcpy(g->pos+3*k,work.pos+3*k,3*sizeof(float)); memcpy(g->nrm+3*k,work.nrm+3*k,3*sizeof(float));
        memcpy(g->tu+3*k,work.tu+3*k,3*sizeof(float)); memcpy(g->tv+3*k,work.tv+3*k,3*sizeof(float));
        g->dist[k] = work.dist[k]; g->support[k] = work.support[k]; g->state[k] = work.state[k];
        g->chart[k] = work.chart[k]; g->component[k] = work.state[k] ? labels[work.component[k]] : -1;
        if (!work.state[k]) out.component_group[k] = -1;
    }
    st->proposed_patches = out.patches; st->proposed_nodes = out.nodes; st->proposed_faces = out.faces;
    fprintf(stderr,"  DATA rebuild: retained %zu DATA; replaced %zu old FILL; %zu nodes/%zu faces in %zu pending patches, "
            "%zu refused (%zu reach, %zu orientation, %zu contact), %zu separated cells, %zu alternate inferred diagonals\n",
            st->data_nodes,st->old_fill_nodes,out.nodes,out.faces,out.patches,st->refused_patches,
            st->reach_refused,st->orientation_refused,st->contact_refused,st->separated_cells,st->alternate_diagonals);
    *result = out; memset(&out,0,sizeof out); rc = 0;
done:
    CloseupGrid_free(&work); CloseupGrid_free(&audit); CloseupRebuildResult_free(&out);
    free(seed); free(mask); free(seen); free(refused); free(d2); free(reference); free(target); free(tu); free(tv); free(weight);
    free(nearest); free(labels); free(parent); free(queue); free(index); free(remap); free(winding);
    return rc;
}

static int cs_rebuild_test_faces(const CloseupGrid *g, int reverse, int32_t **faces, size_t *nf)
{
    size_t nn = g->ni*g->nj;
    *faces = (int32_t *)malloc(6*nn*sizeof(int32_t)); *nf = 0;
    if (!*faces) return -1;
    for (size_t j = 0; j+1 < g->nj; j++) for (size_t i = 0; i+1 < g->ni; i++) {
        size_t a = j*g->ni+i, corners[4] = {a,a+1,a+g->ni,a+g->ni+1};
        /* Alternate the diagonal to exercise retained source topology. */
        int triangles[2] = {2,1};
        if ((i+j)%2) { triangles[0] = 3; triangles[1] = 0; }
        for (int t = 0; t < 2; t++) {
            size_t v[3]; int valid = 1;
            for (int k = 0; k < 3; k++) {
                v[k] = corners[cs_rebuild_triangles[triangles[t]][k]];
                if (g->state[v[k]] != CLOSEUP_NODE_DATA || g->component[v[k]] != g->component[v[0]]) valid = 0;
            }
            if (!valid) continue;
            (*faces)[3*(*nf)] = (int32_t)v[0];
            (*faces)[3*(*nf)+1] = (int32_t)v[reverse ? 2 : 1];
            (*faces)[3*(*nf)+2] = (int32_t)v[reverse ? 1 : 2]; (*nf)++;
        }
    }
    return 0;
}

static int cs_selftest_rebuild(void)
{
    int fails = cs_rebuild_transport_selftest()+cs_rebuild_regrowth_domain_selftest();
    {
        CloseupGrid g = {0}; CloseupRebuildResult proposal = {0}; CloseupRebuildStats st = {0};
        int32_t nearest[9] = {0,1,2,3,0,5,6,7,8}, patch[9] = {-1,-1,-1,-1,0,-1,-1,-1,-1};
        int32_t faces[24] = {0,1,4, 0,4,3, 1,2,5, 1,5,4, 3,4,7, 3,7,6, 4,5,8, 4,8,7};
        float target[27], reference[27] = {0};
        if (CloseupGrid_alloc(&g,-2,2,-2,2,2)) return 1;
        for (size_t k = 0; k < 9; k++) {
            g.state[k] = k == 4 ? CLOSEUP_NODE_FILL : CLOSEUP_NODE_DATA;
            g.pos[3*k] = (float)(12000+g.u0+g.h*(double)(k%g.ni));
            g.pos[3*k+1] = (float)(8000+g.v0+g.h*(double)(k/g.ni)); g.pos[3*k+2] = 4000;
            reference[3*k+2] = 1;
        }
        memcpy(target,g.pos,sizeof target); g.pos[12] += 6;
        proposal.faces = 8; proposal.face_nodes = faces; proposal.node_patch = patch;
        if (cs_rebuild_motion(&g,&proposal,target,nearest,reference,&st) || !st.orientation_before_motion ||
            st.orientation_after_motion || !st.motion_vertices) fails++;
        for (size_t k = 0; k < 9; k++) if (k != 4 && memcmp(g.pos+3*k,target+3*k,3*sizeof(float))) fails++;
        if (fabs((double)g.pos[12]-target[12])>.01 || fabs((double)g.pos[13]-target[13])>.01) fails++;
        fprintf(stderr,"  DATA rebuild inverted-fan constrained motion / exact DATA: %s\n",fails ? "FAIL" : "ok");
        CloseupGrid_free(&g);
    }
    {
        CloseupGrid g = {0}; CloseupRebuildStats st = {0};
        int32_t nearest[4] = {0,1,2,3}; float reference[12] = {0}; int triangles[2] = {0};
        const float point[12] = {0,0,0, 2,0,0, 1.6f,.8f,0, 2,2,0};
        if (CloseupGrid_alloc(&g,0,2,0,2,2)) return 1;
        memcpy(g.pos,point,sizeof point);
        for (size_t k = 0; k < 4; k++) { g.state[k] = CLOSEUP_NODE_DATA; g.component[k] = 0; reference[3*k+2] = 1; }
        if (cs_rebuild_choose_diagonal(&g,0,0,nearest,reference,triangles,&st) != 2 ||
            triangles[0] != 3 || triangles[1] != 0 || st.alternate_diagonals != 1) fails++;
        /* A retained ABD source triangle forbids switching its diagonal, even
         * though its inferred ADC complement will correctly fail orientation. */
        size_t adc[3] = {0,3,2};
        if (cs_rebuild_choose_diagonal(&g,0,4,nearest,reference,triangles,&st) != 2 ||
            triangles[0] != 2 || triangles[1] != 1 || st.alternate_diagonals != 1 ||
            cs_rebuild_triangle_oriented(&g,adc,nearest,reference,NULL)) fails++;
        fprintf(stderr,"  DATA rebuild inferred diagonal / observed diagonal preservation: %s\n",fails ? "FAIL" : "ok");
        CloseupGrid_free(&g);
    }
    for (int test = 0; test < 12; test++) {
        CloseupGrid g = {0}, original = {0}; CloseupRebuildResult out = {0}; CloseupRebuildStats st = {0};
        uint8_t *domain = NULL; int32_t *faces = NULL; size_t nf = 0;
        double extent = test == 6 ? 270 : test == 5 ? 24 : 12;
        int before = fails;
        if (CloseupGrid_alloc(&g,-extent,test == 9 ? 72 : extent,-extent,extent,2)) return 1;
        size_t nn = g.ni*g.nj;
        domain = (uint8_t *)calloc(nn,1);
        if (!domain) { CloseupGrid_free(&g); return 1; }
        for (size_t k = 0; k < nn; k++) {
            double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
            int remote = test == 9 && u >= 58 && u <= 66 && fabs(v)<=4;
            if (!remote && (fabs(u)>extent-2 || fabs(v)>extent-2)) continue;
            if (test == 10 && !(u<=-6 || (u>=2 && u<=6 && fabs(v)<=4))) continue;
            domain[k] = 1;
            int observed = fabs(u)>=extent-4 || fabs(v)>=extent-4;
            if (test < 2) observed = u <= -4 || v <= -4;
            if (test == 2 || test == 3 || test == 4) observed = fabs(u)>=4;
            if (test == 5) observed = ((k%g.ni)%5)<2 || fabs(u)>=extent-4;
            if (test == 7) observed = 0;
            if (test == 10) observed = u<=-6;
            if (test == 11) observed = fabs(u)>=1 || fabs(v)>=1;
            if (remote) observed = 1;
            g.state[k] = observed ? CLOSEUP_NODE_DATA : CLOSEUP_NODE_FILL;
            g.component[k] = ((test == 2 || test == 3) && u>=4) || (test == 10 && !observed) ? 29 : remote ? 41 : 17;
            g.chart[k] = observed ? remote ? 23 : 19 : -1;
            g.support[k] = observed ? 12 : 0; g.dist[k] = observed ? 0 : 6;
            double x = test == 5 ? 10*sin(u/10) : remote ? u-62 : u;
            double z = test == 5 ? 10*cos(u/10) : test == 3 && u>=4 ? 6 : remote ? .5 : 0;
            if (test == 4 && u>0) x -= 16;
            g.pos[3*k] = (float)(12000+x); g.pos[3*k+1] = (float)(8000+v); g.pos[3*k+2] = (float)(4000+z);
            g.tu[3*k] = test == 5 ? (float)cos(u/10) : 1;
            g.tu[3*k+2] = test == 5 ? (float)-sin(u/10) : 0; g.tv[3*k+1] = 1;
            cs_set_normal(&g,k);
            if (test == 1) for (int d = 0; d < 3; d++) g.nrm[3*k+(size_t)d] = -g.nrm[3*k+(size_t)d];
            if (!observed) {
                if (test == 11) g.pos[3*k] += 6; /* an inverted old generated fan */
                else g.pos[3*k+2] += 30; /* erroneous prior inference */
            }
        }
        if (cs_rebuild_test_faces(&g,test == 1,&faces,&nf) || cs_rebuild_copy_grid(&g,&original)) { fails++; goto test_done; }
        if (CloseupSurface_rebuild_generated(&g,domain,faces,nf,test == 8 ? 1 : 1000,&out,&st)) { fails++; goto test_done; }
        if (!st.metric_qualification_pending) fails++;
        size_t generated = 0;
        for (size_t k = 0; k < nn; k++) {
            if (original.state[k] == CLOSEUP_NODE_DATA) {
                if (g.state[k] != original.state[k] || g.chart[k] != original.chart[k] || g.component[k] != original.component[k] ||
                    memcmp(g.pos+3*k,original.pos+3*k,3*sizeof(float)) || memcmp(g.nrm+3*k,original.nrm+3*k,3*sizeof(float)) ||
                    memcmp(g.tu+3*k,original.tu+3*k,3*sizeof(float)) || memcmp(g.tv+3*k,original.tv+3*k,3*sizeof(float)) ||
                    memcmp(g.dist+k,original.dist+k,sizeof(float)) || memcmp(g.support+k,original.support+k,sizeof(float))) fails++;
            } else if (g.state[k]) {
                generated++;
                if (!domain[k] || g.state[k] != CLOSEUP_NODE_FILL || g.chart[k] != -1 || out.node_patch[k] < 0) fails++;
                if (test < 3 || test == 6 || test == 10 || test == 11) {
                    double u = g.u0+g.h*(double)(k%g.ni), v = g.v0+g.h*(double)(k/g.ni);
                    if (fabs((double)g.pos[3*k]-(12000+u))>.01 || fabs((double)g.pos[3*k+1]-(8000+v))>.01 ||
                        fabs((double)g.pos[3*k+2]-4000)>.01 || (test == 1 ? g.nrm[3*k+2]>=0 : g.nrm[3*k+2]<=0)) fails++;
                }
            }
        }
        for (size_t f = 0; f < out.faces; f++) {
            int32_t a = out.face_nodes[3*f], b = out.face_nodes[3*f+1], c = out.face_nodes[3*f+2];
            double det = ((double)((size_t)b%g.ni)-(double)((size_t)a%g.ni))*((double)((size_t)c/g.ni)-(double)((size_t)a/g.ni))-
                         ((double)((size_t)b/g.ni)-(double)((size_t)a/g.ni))*((double)((size_t)c%g.ni)-(double)((size_t)a%g.ni));
            if (test == 1 ? !(det < 0) : !(det > 0)) fails++;
            if (out.component_group[a] != out.component_group[b] || out.component_group[a] != out.component_group[c]) fails++;
        }
        if (test < 3 || test == 5 || test == 6 || test == 11) {
            if (generated != st.old_fill_nodes || st.refused_patches || !out.faces) fails++;
        }
        if (test == 10 && (generated || !st.transport_unanchored_faces || !st.orientation_refused)) fails++;
        if (test == 2 && st.certified_joins != 1) fails++;
        if (test == 3 && (st.certified_joins || st.refused_joins != 1 || !st.separated_cells)) fails++;
        if (test == 4 && !st.orientation_refused && !st.contact_refused) fails++;
        if (test == 6 && generated <= CS_INFILL_MAX_NODES) fails++;
        if ((test == 7 || test == 8) && (generated || !st.refused_patches)) fails++;
        if (test == 9 && !st.contact_refused) fails++;
        if (test == 11 && (!st.target_start_patches || st.old_fill_start_nodes || st.orientation_after_motion)) fails++;
test_done:
        fprintf(stderr,"  DATA rebuild case %d: %zu retained, %zu rebuilt, %zu new faces; %s\n",
                test,st.data_nodes,out.nodes,out.faces,fails == before ? "ok" : "FAIL");
        free(domain); free(faces); CloseupRebuildResult_free(&out); CloseupGrid_free(&g); CloseupGrid_free(&original);
    }
    fprintf(stderr,"  DATA-preserving rebuild: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails != 0;
}

int CloseupSurface_rebuild_selftest(void)
{
    return cs_selftest_rebuild();
}
