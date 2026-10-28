/* Local orientation is a continuation property, not agreement with a normal
 * at a possibly distant DATA donor. An inferred face must belong to an actual
 * edge-connected component containing retained observed faces. A negative
 * normal dot across an edge is a possible sharp crease, not an inversion.
 * Complete existing 3-D intersection predicates distinguish true fold-backs
 * and crossings from those creases. Final UV/material/metric/signal checks
 * remain separate, mandatory gates. */
#include "../remesh/intersection_cleanup.h"

typedef struct CsRebuildTransport {
    size_t components, anchored_components, unanchored_components;
    size_t unanchored_faces, degenerate_faces, bad_edges, winding_edges, nonmanifold_edges;
    CloseupRebuildWorldStats exact, proximity;
} CsRebuildTransport;

typedef struct CsRebuildWorldVisit {
    const CloseupRebuildResult *proposal;
    size_t observed_nf;
    uint8_t *refused;
    CloseupRebuildWorldStats *stats;
} CsRebuildWorldVisit;

static int cs_rebuild_world_visit(size_t a, size_t b, int kind, void *context)
{
    CsRebuildWorldVisit *v = (CsRebuildWorldVisit *)context;
    v->stats->pairs++;
    v->stats->overlaps += kind == INTERSECTION_HIT_OVERLAP;
    v->stats->stabs += kind == INTERSECTION_HIT_STAB;
    v->stats->folds += kind == INTERSECTION_HIT_FOLD;
    if (a < v->observed_nf && b < v->observed_nf) v->stats->observed_observed++;
    else if (a >= v->observed_nf && b >= v->observed_nf) v->stats->inferred_inferred++;
    else v->stats->observed_inferred++;
    if (a >= v->observed_nf) {
        int32_t patch = v->proposal->face_patch[a-v->observed_nf];
        if (!v->refused[patch]) v->refused[patch] = 3;
    }
    if (b >= v->observed_nf) {
        int32_t patch = v->proposal->face_patch[b-v->observed_nf];
        if (!v->refused[patch]) v->refused[patch] = 3;
    }
    return 0;
}

static const int32_t *cs_rebuild_transport_face(const CloseupRebuildResult *proposal,
                                               const int32_t *observed, size_t observed_nf, size_t f)
{
    return f < observed_nf ? observed+3*f : proposal->face_nodes+3*(f-observed_nf);
}

static void cs_rebuild_transport_refuse(const CloseupRebuildResult *proposal, size_t observed_nf,
                                        size_t f, uint8_t *refused)
{
    if (f >= observed_nf) refused[proposal->face_patch[f-observed_nf]] = 1;
}

static int cs_rebuild_transport_certificate(const CloseupGrid *g, const CloseupRebuildResult *proposal,
                                             const int32_t *observed, size_t observed_nf,
                                             uint8_t *refused, CsRebuildTransport *st)
{
    size_t nn = g->ni*g->nj, nf = observed_nf+proposal->faces;
    int32_t *edge = NULL, *parent = NULL, *all_faces = NULL;
    int8_t *direction = NULL;
    uint8_t *multiplicity = NULL, *anchored = NULL;
    double *normal = NULL;
    int rc = -1;
    memset(st,0,sizeof *st);
    if (nf > INT32_MAX || nn > SIZE_MAX/(4*sizeof(int32_t))) return -1;
    if (!nf) { st->exact.complete = st->proximity.complete = 1; return 0; }
    edge = (int32_t *)malloc(4*nn*sizeof(int32_t));
    direction = (int8_t *)calloc(4*nn,1); multiplicity = (uint8_t *)calloc(4*nn,1);
    parent = (int32_t *)malloc((nf ? nf : 1)*sizeof(int32_t));
    anchored = (uint8_t *)calloc(nf ? nf : 1,1);
    normal = (double *)calloc(nf ? 3*nf : 1,sizeof(double));
    all_faces = (int32_t *)malloc(3*nf*sizeof(int32_t));
    if (!edge || !direction || !multiplicity || !parent || !anchored || !normal || !all_faces) goto done;
    if (observed_nf) memcpy(all_faces,observed,3*observed_nf*sizeof(int32_t));
    if (proposal->faces) memcpy(all_faces+3*observed_nf,proposal->face_nodes,3*proposal->faces*sizeof(int32_t));
    for (size_t k = 0; k < 4*nn; k++) edge[k] = -1;
    for (size_t f = 0; f < nf; f++) {
        const int32_t *v = cs_rebuild_transport_face(proposal,observed,observed_nf,f);
        double a[3], b[3], length2 = 0;
        parent[f] = (int32_t)f;
        for (int k = 0; k < 3; k++) if (v[k] < 0 || (size_t)v[k] >= nn) goto done;
        for (int d = 0; d < 3; d++) {
            a[d] = (double)g->pos[3*(size_t)v[1]+(size_t)d]-g->pos[3*(size_t)v[0]+(size_t)d];
            b[d] = (double)g->pos[3*(size_t)v[2]+(size_t)d]-g->pos[3*(size_t)v[0]+(size_t)d];
        }
        for (int d = 0; d < 3; d++) {
            double value = a[(d+1)%3]*b[(d+2)%3]-a[(d+2)%3]*b[(d+1)%3];
            normal[3*f+(size_t)d] = value; length2 += value*value;
        }
        if (!(length2 > 0) || !isfinite(length2)) {
            if (f < observed_nf) goto done;
            st->degenerate_faces++; cs_rebuild_transport_refuse(proposal,observed_nf,f,refused);
        } else for (int d = 0; d < 3; d++) normal[3*f+(size_t)d] /= sqrt(length2);
    }
    for (size_t f = 0; f < nf; f++) {
        const int32_t *v = cs_rebuild_transport_face(proposal,observed,observed_nf,f);
        for (int k = 0; k < 3; k++) {
            size_t a = (size_t)v[k], b = (size_t)v[(k+1)%3]; int8_t sign = 1;
            if (a > b) { size_t swap = a; a = b; b = swap; sign = -1; }
            size_t dx = b%g->ni > a%g->ni ? b%g->ni-a%g->ni : a%g->ni-b%g->ni;
            size_t dy = b/g->ni-a/g->ni, delta = b-a;
            if (dx > 1 || dy > 1 || (!dx && !dy)) goto done;
            size_t kind = delta == 1 ? 0 : delta == g->ni ? 1 : delta == g->ni+1 ? 2 : 3;
            size_t key = 4*a+kind;
            if (edge[key] < 0) { edge[key] = (int32_t)f; direction[key] = sign; multiplicity[key] = 1; continue; }
            size_t other = (size_t)edge[key];
            if (multiplicity[key] >= 2) {
                if (multiplicity[key] == 2) st->nonmanifold_edges++;
                multiplicity[key] = 3;
                cs_rebuild_transport_refuse(proposal,observed_nf,f,refused);
                cs_rebuild_transport_refuse(proposal,observed_nf,other,refused);
                continue;
            }
            multiplicity[key] = 2;
            int32_t first = cs_rebuild_find(parent,(int32_t)f), second = cs_rebuild_find(parent,(int32_t)other);
            if (first != second) parent[second] = first;
            if (direction[key] == sign) {
                st->winding_edges++;
                cs_rebuild_transport_refuse(proposal,observed_nf,f,refused);
                cs_rebuild_transport_refuse(proposal,observed_nf,other,refused);
            }
            if (f >= observed_nf || other >= observed_nf) {
                double cosine = 0;
                for (int d = 0; d < 3; d++) cosine += normal[3*f+(size_t)d]*normal[3*other+(size_t)d];
                if (!(cosine > 0) || !isfinite(cosine)) {
                    st->bad_edges++;
                }
            }
        }
    }
    for (size_t f = 0; f < observed_nf; f++) anchored[cs_rebuild_find(parent,(int32_t)f)] = 1;
    for (size_t f = 0; f < nf; f++) if (parent[f] == (int32_t)f) {
        st->components++; st->anchored_components += anchored[f] != 0; st->unanchored_components += !anchored[f];
    }
    for (size_t f = observed_nf; f < nf; f++) if (!anchored[cs_rebuild_find(parent,(int32_t)f)]) {
        st->unanchored_faces++; cs_rebuild_transport_refuse(proposal,observed_nf,f,refused);
    }
    /* Never pass an inferred-only mask: that API requires BOTH faces marked
     * and would omit every observed/inferred collision. Existing observed
     * conflicts are counted separately and never deleted or hidden. */
    for (int mode = 0; mode < 2; mode++) {
        IntersectionCleanupParams params; IntersectionCleanupStats measured = {0};
        CloseupRebuildWorldStats *stats = mode ? &st->proximity : &st->exact;
        CsRebuildWorldVisit visitor = {proposal,observed_nf,refused,stats};
        IntersectionCleanup_default_params(&params);
        if (!mode) params.gap_max = 0;
        params.max_conflicts = SIZE_MAX;
        if (IntersectionCleanup_audit_visit_parallel(g->pos,nn,all_faces,nf,NULL,&params,NULL,
                                                     cs_rebuild_world_visit,&visitor,&measured)) goto done;
        if (stats->pairs != measured.conflicts) goto done;
        stats->complete = 1;
    }
    rc = 0;
done:
    free(edge); free(direction); free(multiplicity); free(parent); free(anchored); free(normal); free(all_faces); return rc;
}

static int cs_rebuild_transport_selftest(void)
{
    int fails = 0;
    /* A one-sided observed collar supports a 229-degree smooth cylinder.
     * The same topology folded back through 180 degrees is locally invalid;
     * a detached smooth inferred island has no observed orientation anchor. */
    for (int test = 0; test < 5; test++) {
        CloseupGrid g = {0}; CloseupRebuildResult proposal = {0}; CsRebuildTransport st = {0};
        int32_t *faces = NULL; uint8_t refused[1] = {0}; size_t nf = 0, observed_nf = 0;
        int before = fails;
        if (CloseupGrid_alloc(&g,0,80,0,4,2)) return 1;
        size_t nn = g.ni*g.nj, maxfaces = 2*(g.ni-1)*(g.nj-1);
        faces = (int32_t *)malloc(3*maxfaces*sizeof(int32_t));
        proposal.face_patch = (int32_t *)calloc(maxfaces,sizeof(int32_t));
        if (!faces || !proposal.face_patch) { free(faces); free(proposal.face_patch); CloseupGrid_free(&g); return 1; }
        for (size_t k = 0; k < nn; k++) {
            double u = g.h*(double)(k%g.ni), v = g.h*(double)(k/g.ni);
            double angle = u/20;
            g.state[k] = k%g.ni < 2 && test != 2 ? CLOSEUP_NODE_DATA : CLOSEUP_NODE_FILL;
            g.pos[3*k] = test == 1 ? (float)(12000+(u <= 40 ? u : 80-u)) : (float)(12000+20*sin(angle));
            g.pos[3*k+1] = (float)(8000+v);
            g.pos[3*k+2] = test == 1 ? 4000 : (float)(4000+20*cos(angle));
            if (test == 4) {
                double theta = 110*acos(-1.0)/180, along = fmax(0,u-40);
                g.pos[3*k] = (float)(12000+fmin(u,40)+along*cos(theta));
                g.pos[3*k+2] = (float)(4000+along*sin(theta));
            }
        }
        /* Observed faces are an exact prefix; both winding choices work. */
        for (int pass = 0; pass < 2; pass++) for (size_t j = 0; j+1 < g.nj; j++) for (size_t i = 0; i+1 < g.ni; i++) {
            int observed = i == 0 && test != 2;
            if ((pass == 0) != observed) continue;
            size_t a = j*g.ni+i; const size_t tri[6] = {a,a+1,a+g.ni+1,a,a+g.ni+1,a+g.ni};
            for (int t = 0; t < 2; t++) {
                faces[3*nf] = (int32_t)tri[3*t];
                faces[3*nf+1] = (int32_t)tri[3*t+(test == 3 ? 2 : 1)];
                faces[3*nf+2] = (int32_t)tri[3*t+(test == 3 ? 1 : 2)]; nf++;
                observed_nf += observed != 0;
            }
        }
        proposal.face_nodes = faces+3*observed_nf; proposal.faces = nf-observed_nf;
        if (cs_rebuild_transport_certificate(&g,&proposal,faces,observed_nf,refused,&st)) fails++;
        if (test == 0 || test == 3 || test == 4) {
            if (refused[0] || st.unanchored_faces || st.anchored_components != 1 || st.exact.pairs || st.proximity.pairs) fails++;
            if (test == 4 && !st.bad_edges) fails++;
        } else if (test == 1) { if (!refused[0] || !st.bad_edges || st.unanchored_faces || !st.exact.pairs) fails++; }
        else if (!refused[0] || st.bad_edges || st.unanchored_faces != proposal.faces) fails++;
        if (!st.exact.complete || !st.proximity.complete) fails++;
        fprintf(stderr,"  DATA rebuild local orientation case %d (cylinder/fold/unanchored/reverse/110-degree crease): %s\n",test,fails == before ? "ok" : "FAIL");
        free(faces); free(proposal.face_patch); CloseupGrid_free(&g);
    }
    {
        const float xyz[18] = {0,0,0, 4,0,0, 0,4,0, 1,1,-2, 1,1,2, 3,1,0};
        const int32_t a[3] = {0,1,2}, b[3] = {3,4,5};
        IntersectionCleanupParams params; IntersectionCleanup_default_params(&params); params.gap_max = 0;
        if (IntersectionCleanup_pair_test(xyz,6,a,b,&params,NULL) != INTERSECTION_HIT_STAB) fails++;
    }
    return fails != 0;
}
