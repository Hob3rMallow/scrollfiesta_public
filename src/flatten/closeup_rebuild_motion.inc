/* Generated-only constraint projections toward the Cartesian TAUCS target.
 * A triangle's signed projected area is affine in one moving vertex. Grid
 * colors ensure simultaneous moves never share a triangle. These are local
 * 3-D feasibility operations, not an alternative UV flattening algorithm. */
typedef struct CsRebuildConstraint {
    int32_t face;
    double area, margin, gradient[3], norm2;
} CsRebuildConstraint;

static int cs_rebuild_choose_start(CloseupGrid *g, const CloseupGrid *original,
                                    const CloseupRebuildResult *proposal, size_t npatches,
                                    const float *target, const int32_t *nearest,
                                    const float *reference, CloseupRebuildStats *st)
{
    size_t *old_bad = (size_t *)calloc(npatches ? npatches : 1,sizeof(size_t));
    size_t *target_bad = (size_t *)calloc(npatches ? npatches : 1,sizeof(size_t));
    CloseupGrid target_grid = *g; target_grid.pos = (float *)target;
    if (!old_bad || !target_bad) { free(old_bad); free(target_bad); return -1; }
    for (size_t f = 0; f < proposal->faces; f++) {
        size_t v[3];
        for (int k = 0; k < 3; k++) v[k] = (size_t)proposal->face_nodes[3*f+(size_t)k];
        double det = ((double)(v[1]%g->ni)-(double)(v[0]%g->ni))*((double)(v[2]/g->ni)-(double)(v[0]/g->ni))-
                     ((double)(v[1]/g->ni)-(double)(v[0]/g->ni))*((double)(v[2]%g->ni)-(double)(v[0]%g->ni));
        if (det < 0) { size_t swap = v[1]; v[1] = v[2]; v[2] = swap; }
        int32_t patch = proposal->face_patch[f];
        old_bad[patch] += !cs_rebuild_triangle_oriented(g,v,nearest,reference,NULL);
        target_bad[patch] += !cs_rebuild_triangle_oriented(&target_grid,v,nearest,reference,NULL);
    }
    st->old_fill_start_nodes = 0;
    for (size_t p = 0; p < npatches; p++) st->target_start_patches += target_bad[p] <= old_bad[p];
    for (size_t k = 0; k < g->ni*g->nj; k++) if (proposal->node_patch[k] >= 0) {
        int32_t p = proposal->node_patch[k];
        if (target_bad[p] <= old_bad[p]) memcpy(g->pos+3*k,target+3*k,3*sizeof(float));
        else st->old_fill_start_nodes += original->state[k] == CLOSEUP_NODE_FILL;
    }
    free(old_bad); free(target_bad); return 0;
}

static double cs_rebuild_motion_area(const CloseupGrid *g, const int32_t *face, const double *normal)
{
    double e[3], f[3], area = 0;
    for (int d = 0; d < 3; d++) {
        e[d] = (double)g->pos[3*(size_t)face[1]+(size_t)d]-g->pos[3*(size_t)face[0]+(size_t)d];
        f[d] = (double)g->pos[3*(size_t)face[2]+(size_t)d]-g->pos[3*(size_t)face[0]+(size_t)d];
    }
    for (int d = 0; d < 3; d++) area += normal[d]*(e[(d+1)%3]*f[(d+2)%3]-e[(d+2)%3]*f[(d+1)%3]);
    return area;
}

static int cs_rebuild_motion_vertex(CloseupGrid *g, size_t at, const float *target,
                                     const int32_t *faces, const double *normal,
                                     const size_t *offset, const int32_t *incident)
{
    CsRebuildConstraint constraints[8];
    size_t count = offset[at+1]-offset[at], bad = 0;
    double before_bad = 0, before_target = 0, edge_scale2 = 0, delta[3] = {0};
    float old[3];
    if (!count || count > 8) return 0;
    memcpy(old,g->pos+3*at,sizeof old);
    for (size_t q = 0; q < count; q++) {
        int32_t f = incident[offset[at]+q]; const int32_t *face = faces+3*(size_t)f;
        const double *n = normal+3*(size_t)f;
        int corner = face[0] == (int32_t)at ? 0 : face[1] == (int32_t)at ? 1 : 2;
        size_t b = (size_t)face[(corner+1)%3], c = (size_t)face[(corner+2)%3];
        double opposite[3], roundoff = 0;
        CsRebuildConstraint *constraint = constraints+q;
        constraint->face = f; constraint->area = cs_rebuild_motion_area(g,face,n); constraint->norm2 = 0;
        for (int d = 0; d < 3; d++) {
            opposite[d] = (double)g->pos[3*b+(size_t)d]-g->pos[3*c+(size_t)d];
            double ab = (double)g->pos[3*b+(size_t)d]-old[d], ac = (double)g->pos[3*c+(size_t)d]-old[d];
            edge_scale2 += ab*ab+ac*ac;
        }
        for (int d = 0; d < 3; d++) {
            double gradient = opposite[(d+1)%3]*n[(d+2)%3]-opposite[(d+2)%3]*n[(d+1)%3];
            double up = fabs((double)nextafterf(old[d],INFINITY)-old[d]);
            double down = fabs((double)nextafterf(old[d],-INFINITY)-old[d]);
            constraint->gradient[d] = gradient; constraint->norm2 += gradient*gradient;
            roundoff += 2*fabs(gradient)*fmax(up,down);
        }
        constraint->margin = constraint->area > 0 ? fmin(.5*constraint->area,roundoff) : roundoff;
        if (!(constraint->area > 0)) {
            bad++;
            double violation = constraint->margin-constraint->area;
            if (constraint->norm2 > 0) before_bad += violation*violation/constraint->norm2;
        }
    }
    for (int d = 0; d < 3; d++) {
        double desired = (double)target[3*at+(size_t)d]-old[d];
        before_target += desired*desired;
        if (!bad) delta[d] = desired;
    }
    /* Alternating half-space projections restore invalid local areas first.
     * Once feasible, they move toward the variational target. No valid face
     * is allowed to become invalid during either phase. */
    for (int sweep = 0; sweep < CLOSEUP_REBUILD_PROJECTION_SWEEPS; sweep++) {
        int projected = 0;
        for (size_t q = 0; q < count; q++) {
            const CsRebuildConstraint *constraint = constraints+q;
            double predicted = constraint->area;
            for (int d = 0; d < 3; d++) predicted += constraint->gradient[d]*delta[d];
            if (predicted >= constraint->margin || !(constraint->norm2 > 0)) continue;
            double step = (constraint->margin-predicted)/constraint->norm2;
            for (int d = 0; d < 3; d++) delta[d] += step*constraint->gradient[d];
            projected = 1;
        }
        if (!projected) break;
    }
    double length2 = 0, alpha = 1;
    for (int d = 0; d < 3; d++) length2 += delta[d]*delta[d];
    /* A fully collapsed triangle has zero single-vertex area gradients. A
     * target-directed neutral step can separate its coincident vertices,
     * after which the next color can restore strictly positive area. */
    if (!(length2 > 0) && bad) for (int d = 0; d < 3; d++) {
        delta[d] = (double)target[3*at+(size_t)d]-old[d]; length2 += delta[d]*delta[d];
    }
    double trust2 = edge_scale2/(2*(double)count);
    if (trust2 > 0 && length2 > trust2) alpha = sqrt(trust2/length2);
    for (size_t q = 0; q < count; q++) if (constraints[q].area > 0) {
        double change = 0;
        for (int d = 0; d < 3; d++) change += constraints[q].gradient[d]*delta[d];
        if (change < 0) alpha = fmin(alpha,.9*constraints[q].area/-change);
    }
    if (!(alpha > 0) || !isfinite(alpha)) return 0;
    for (int d = 0; d < 3; d++) g->pos[3*at+(size_t)d] = (float)((double)old[d]+alpha*delta[d]);
    if (!memcmp(g->pos+3*at,old,sizeof old)) return 0;
    size_t after_bad_count = 0; double after_bad = 0, after_target = 0;
    int preserve = 1;
    for (size_t q = 0; q < count; q++) {
        const CsRebuildConstraint *constraint = constraints+q;
        size_t f = (size_t)constraint->face;
        double area = cs_rebuild_motion_area(g,faces+3*f,normal+3*f);
        if (!isfinite(area) || (constraint->area > 0 && !(area > 0))) preserve = 0;
        if (!(area > 0)) {
            double violation = fmax(0,constraint->margin-area);
            after_bad_count++;
            if (constraint->norm2 > 0) after_bad += violation*violation/constraint->norm2;
        }
    }
    for (int d = 0; d < 3; d++) {
        double error = (double)target[3*at+(size_t)d]-g->pos[3*at+(size_t)d];
        after_target += error*error;
        if (!isfinite(g->pos[3*at+(size_t)d])) preserve = 0;
    }
    int progress = bad ? after_bad_count < bad ||
        (after_bad_count == bad && (after_bad < before_bad || (after_bad <= before_bad && after_target < before_target))) :
        after_target < before_target;
    if (!preserve || !progress) { memcpy(g->pos+3*at,old,sizeof old); return 0; }
    return 1;
}

static size_t cs_rebuild_motion_measure(const CloseupGrid *g, const CloseupRebuildResult *proposal,
                                        const float *target, const double *normal,
                                        double *bad_energy, double *target_error)
{
    size_t bad = 0; *bad_energy = *target_error = 0;
    for (size_t f = 0; f < proposal->faces; f++) {
        double area = cs_rebuild_motion_area(g,proposal->face_nodes+3*f,normal+3*f);
        if (!(area > 0)) { bad++; *bad_energy += area*area; }
    }
    for (size_t k = 0; k < g->ni*g->nj; k++) if (proposal->node_patch[k] >= 0) for (int d = 0; d < 3; d++) {
        double error = (double)g->pos[3*k+(size_t)d]-target[3*k+(size_t)d]; *target_error += error*error;
    }
    return bad;
}

#include "closeup_rebuild_coupled.inc"

static int cs_rebuild_motion(CloseupGrid *g, const CloseupRebuildResult *proposal, const float *target,
                              const int32_t *nearest, const float *reference, CloseupRebuildStats *st)
{
    size_t nn = g->ni*g->nj;
    size_t *offset = (size_t *)calloc(nn+1,sizeof(size_t)), *cursor = (size_t *)malloc(nn*sizeof(size_t));
    int32_t *incident = (int32_t *)malloc((proposal->faces ? 3*proposal->faces : 1)*sizeof(int32_t));
    double *normal = (double *)calloc(proposal->faces ? 3*proposal->faces : 1,sizeof(double));
    CsRebuildBlock *block = (CsRebuildBlock *)malloc(sizeof *block);
    int rc = -1;
    if (!offset || !cursor || !incident || !normal || !block) goto done;
    for (size_t f = 0; f < proposal->faces; f++) {
        const int32_t *face = proposal->face_nodes+3*f;
        double length2 = 0;
        double au = (double)((size_t)face[1]%g->ni)-(double)((size_t)face[0]%g->ni);
        double av = (double)((size_t)face[1]/g->ni)-(double)((size_t)face[0]/g->ni);
        double bu = (double)((size_t)face[2]%g->ni)-(double)((size_t)face[0]%g->ni);
        double bv = (double)((size_t)face[2]/g->ni)-(double)((size_t)face[0]/g->ni);
        double sign = au*bv-av*bu > 0 ? 1 : -1;
        for (int k = 0; k < 3; k++) {
            size_t at = (size_t)face[k]; int32_t data = nearest[at]; offset[at+1]++;
            if (data >= 0) for (int d = 0; d < 3; d++) normal[3*f+(size_t)d] += reference[3*(size_t)data+(size_t)d];
        }
        for (int d = 0; d < 3; d++) length2 += normal[3*f+(size_t)d]*normal[3*f+(size_t)d];
        if (length2 > 0) for (int d = 0; d < 3; d++) normal[3*f+(size_t)d] *= sign/sqrt(length2);
    }
    for (size_t k = 0; k < nn; k++) { offset[k+1] += offset[k]; cursor[k] = offset[k]; }
    for (size_t f = 0; f < proposal->faces; f++) for (int k = 0; k < 3; k++) incident[cursor[proposal->face_nodes[3*f+(size_t)k]]++] = (int32_t)f;
    double bad_energy = 0, target_error = 0;
    size_t bad = cs_rebuild_motion_measure(g,proposal,target,normal,&bad_energy,&target_error);
    st->orientation_before_motion = bad; st->target_error_before_motion = target_error;
    fprintf(stderr,"  DATA rebuild local motion: %zu initial invalid triangles, target error %.9g, %zu old FILL starting nodes\n",
            bad,target_error,st->old_fill_start_nodes);
    for (int sweep = 0; sweep < CLOSEUP_REBUILD_MOTION_SWEEPS; sweep++) {
        size_t moved = 0, prior_bad = bad; double prior_energy = bad_energy, prior_error = target_error;
        for (int color_index = 0; color_index < 4; color_index++) {
            int color = sweep%2 ? 3-color_index : color_index;
            int row = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:moved)
#endif
            for (row = 0; row < (int)g->nj; row++) {
                size_t j = (size_t)row;
                if ((int)(j%2) != color/2) continue;
                for (size_t i = (size_t)(color%2); i < g->ni; i += 2) {
                    size_t at = j*g->ni+i;
                    if (proposal->node_patch[at] < 0) continue;
                    moved += (size_t)cs_rebuild_motion_vertex(g,at,target,proposal->face_nodes,normal,offset,incident);
                }
            }
        }
        bad = cs_rebuild_motion_measure(g,proposal,target,normal,&bad_energy,&target_error);
        st->motion_sweeps++; st->motion_vertices += moved;
        fprintf(stderr,"  DATA rebuild local motion sweep %d: %zu moved, invalid %zu -> %zu, negative-area energy %.9g, target error %.9g\n",
                sweep+1,moved,prior_bad,bad,bad_energy,target_error);
        if (bad && (bad == prior_bad && bad_energy >= prior_energy*(1-FLT_EPSILON))) {
            size_t before_blocks = st->coupled_blocks, before_bad = bad;
            double before_energy = bad_energy;
            st->coupled_passes++;
            for (size_t f = 0; f < proposal->faces; f++) {
                int accepted = cs_rebuild_motion_block(g,proposal,f,target,normal,offset,incident,block,st);
                if (accepted < 0) goto done;
                moved += accepted != 0;
            }
            bad = cs_rebuild_motion_measure(g,proposal,target,normal,&bad_energy,&target_error);
            fprintf(stderr,"  DATA rebuild coupled pass %zu: %zu accepted blocks, invalid %zu -> %zu, negative-area energy %.9g -> %.9g\n",
                    st->coupled_passes,st->coupled_blocks-before_blocks,before_bad,bad,before_energy,bad_energy);
        }
        if (!moved || (bad == prior_bad && bad_energy >= prior_energy*(1-FLT_EPSILON) &&
                      target_error >= prior_error*(1-FLT_EPSILON))) break;
    }
    st->orientation_after_motion = bad; st->target_error_after_motion = target_error;
    rc = 0;
done:
    free(offset); free(cursor); free(incident); free(normal); free(block); return rc;
}
