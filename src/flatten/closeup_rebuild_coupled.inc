/* A triangle and its incident vertex stars fit within a 4-by-4 grid block.
 * Moving those generated nodes together resolves constraints that have no
 * feasible single-vertex update. DATA never enters the variable set. */
enum { CS_REBUILD_BLOCK_NODES = 16, CS_REBUILD_BLOCK_VARS = 48, CS_REBUILD_BLOCK_FACES = 128 };

typedef struct CsRebuildBlock {
    int32_t nodes[CS_REBUILD_BLOCK_NODES], faces[CS_REBUILD_BLOCK_FACES];
    float saved[CS_REBUILD_BLOCK_VARS];
    double gradient[CS_REBUILD_BLOCK_FACES*CS_REBUILD_BLOCK_VARS];
    double area[CS_REBUILD_BLOCK_FACES], margin[CS_REBUILD_BLOCK_FACES];
    double rhs[CS_REBUILD_BLOCK_FACES], multiplier[CS_REBUILD_BLOCK_FACES];
    double delta[CS_REBUILD_BLOCK_VARS], trust2[CS_REBUILD_BLOCK_NODES], norm[CS_REBUILD_BLOCK_FACES];
    int active_row[CS_REBUILD_BLOCK_FACES], column[CS_REBUILD_BLOCK_FACES+1];
    int row[CS_REBUILD_BLOCK_FACES*(CS_REBUILD_BLOCK_FACES+1)/2];
    double value[CS_REBUILD_BLOCK_FACES*(CS_REBUILD_BLOCK_FACES+1)/2];
    uint8_t active[CS_REBUILD_BLOCK_FACES];
} CsRebuildBlock;

static int cs_rebuild_block_index(const int32_t *values, int count, int32_t wanted)
{
    for (int k = 0; k < count; k++) if (values[k] == wanted) return k;
    return -1;
}

/* Earliest strictly positive zero of a + b*t + c*t*t. The starting
 * triangle has a>0; clipping before this zero preserves the entire path. */
static double cs_rebuild_area_first_zero(double a, double b, double c)
{
    if (c == 0) return b < 0 ? -a/b : INFINITY;
    double discriminant = b*b-4*a*c;
    if (!(discriminant >= 0) || !isfinite(discriminant)) return INFINITY;
    double q = -.5*(b+copysign(sqrt(discriminant),b));
    double first = INFINITY, root = q/c;
    if (root > 0) first = root;
    if (q != 0) { root = a/q; if (root > 0 && root < first) first = root; }
    return first;
}

static int cs_rebuild_motion_block(CloseupGrid *g, const CloseupRebuildResult *proposal, size_t seed_face,
                                    const float *target, const double *normal, const size_t *offset,
                                    const int32_t *incident, CsRebuildBlock *block, CloseupRebuildStats *st)
{
    int nodes = 0, faces = 0, variables = 0, rc = 0;
    size_t bad = 0;
    double energy = 0, alpha = 1;
    SparseFactor_T factor = NULL;
    int32_t patch = proposal->face_patch[seed_face];
    const int32_t *seed = proposal->face_nodes+3*seed_face;
    if (cs_rebuild_motion_area(g,seed,normal+3*seed_face) > 0) return 0;
    memset(block,0,sizeof *block);
    for (int corner = 0; corner < 3; corner++) {
        size_t at = (size_t)seed[corner];
        for (size_t q = offset[at]; q < offset[at+1]; q++) {
            const int32_t *face = proposal->face_nodes+3*(size_t)incident[q];
            for (int k = 0; k < 3; k++) {
                int32_t v = face[k];
                if (proposal->node_patch[v] != patch || cs_rebuild_block_index(block->nodes,nodes,v) >= 0) continue;
                if (nodes == CS_REBUILD_BLOCK_NODES) return -1;
                block->nodes[nodes++] = v;
            }
        }
    }
    if (!nodes) return 0;
    variables = 3*nodes;
    for (int k = 0; k < nodes; k++) {
        size_t at = (size_t)block->nodes[k];
        if (g->state[at] != CLOSEUP_NODE_FILL) return -1;
        memcpy(block->saved+3*k,g->pos+3*at,3*sizeof(float));
        for (size_t q = offset[at]; q < offset[at+1]; q++) {
            int32_t f = incident[q]; const int32_t *face = proposal->face_nodes+3*(size_t)f;
            if (cs_rebuild_block_index(block->faces,faces,f) < 0) {
                if (faces == CS_REBUILD_BLOCK_FACES) return -1;
                block->faces[faces++] = f;
            }
            for (int b = 0; b < 3; b++) if (face[b] != (int32_t)at) for (int d = 0; d < 3; d++) {
                double edge = (double)g->pos[3*(size_t)face[b]+(size_t)d]-g->pos[3*at+(size_t)d];
                block->trust2[k] += edge*edge;
            }
        }
        if (offset[at+1] > offset[at]) block->trust2[k] /= 2*(double)(offset[at+1]-offset[at]);
    }
    for (int r = 0; r < faces; r++) {
        size_t f = (size_t)block->faces[r]; const int32_t *face = proposal->face_nodes+3*f;
        const double *n = normal+3*f; double roundoff = 0;
        block->area[r] = cs_rebuild_motion_area(g,face,n);
        block->active[r] = !(block->area[r] > 0);
        if (block->active[r]) { bad++; energy += block->area[r]*block->area[r]; }
        for (int corner = 0; corner < 3; corner++) {
            int v = cs_rebuild_block_index(block->nodes,nodes,face[corner]);
            if (v < 0) continue;
            size_t a = (size_t)face[corner], b = (size_t)face[(corner+1)%3], c = (size_t)face[(corner+2)%3];
            double opposite[3];
            for (int d = 0; d < 3; d++) opposite[d] = (double)g->pos[3*b+(size_t)d]-g->pos[3*c+(size_t)d];
            for (int d = 0; d < 3; d++) {
                double gradient = opposite[(d+1)%3]*n[(d+2)%3]-opposite[(d+2)%3]*n[(d+1)%3];
                float old = g->pos[3*a+(size_t)d];
                double ulp = fmax(fabs((double)nextafterf(old,INFINITY)-old),fabs((double)nextafterf(old,-INFINITY)-old));
                block->gradient[(size_t)r*CS_REBUILD_BLOCK_VARS+3*(size_t)v+(size_t)d] = gradient;
                roundoff += 2*fabs(gradient)*ulp;
            }
        }
        block->margin[r] = block->area[r] > 0 ? fmin(.5*block->area[r],roundoff) : roundoff;
    }
    if (!bad) return 0;
    /* Minimum-motion projection onto the linearized area half-spaces. The
     * working set solves G*G' multipliers through TAUCS; negative multipliers
     * release constraints. Normalizing rows changes no half-space. A small
     * numerical diagonal protects dependent local rows; actual geometry is
     * accepted only by the nonlinear path and float checks below. */
    memset(block->active,0,sizeof block->active);
    for (int r = 0; r < faces; r++) {
        double norm2 = 0;
        for (int v = 0; v < variables; v++) {
            double value = block->gradient[(size_t)r*CS_REBUILD_BLOCK_VARS+(size_t)v]; norm2 += value*value;
        }
        block->norm[r] = sqrt(norm2);
    }
    for (int iteration = 0; iteration < CLOSEUP_REBUILD_PROJECTION_SWEEPS; iteration++) {
        int worst = -1; double violation = 0;
        for (int r = 0; r < faces; r++) if (block->norm[r] > 0) {
            double remaining = block->margin[r]-block->area[r];
            for (int v = 0; v < variables; v++) remaining -= block->gradient[(size_t)r*CS_REBUILD_BLOCK_VARS+(size_t)v]*block->delta[v];
            remaining /= block->norm[r];
            if (remaining > violation) { violation = remaining; worst = r; }
        }
        if (worst < 0 || block->active[worst]) break;
        block->active[worst] = 1;
        for (int release = 0; release < CLOSEUP_REBUILD_PROJECTION_SWEEPS; release++) {
            int active = 0, nz = 0, negative = -1; double most_negative = 0;
            for (int r = 0; r < faces; r++) if (block->active[r]) block->active_row[active++] = r;
            memset(block->delta,0,sizeof block->delta);
            if (!active) break;
            for (int col = 0; col < active; col++) {
                int a = block->active_row[col]; block->column[col] = nz;
                block->rhs[col] = (block->margin[a]-block->area[a])/block->norm[a];
                for (int row = col; row < active; row++) {
                    int b = block->active_row[row]; double value = 0;
                    for (int v = 0; v < variables; v++) value += block->gradient[(size_t)a*CS_REBUILD_BLOCK_VARS+(size_t)v]*
                        block->gradient[(size_t)b*CS_REBUILD_BLOCK_VARS+(size_t)v];
                    value /= block->norm[a]*block->norm[b];
                    if (a == b) value += sqrt(DBL_EPSILON);
                    block->row[nz] = row; block->value[nz++] = value;
                }
            }
            block->column[active] = nz;
            int zero = 0;
            int failed = active == 1 ? Sparse_factor_spd(1,1,&zero,&zero,block->value,&factor) :
                Sparse_factor_spd_ccs(active,block->column,block->row,block->value,&factor);
            if (failed || Sparse_factor_solve(factor,block->rhs,block->multiplier)) { rc = -1; goto done; }
            Sparse_factor_free(&factor); st->coupled_solves++;
            for (int r = 0; r < active; r++) {
                int a = block->active_row[r]; double multiplier = block->multiplier[r];
                if (multiplier < most_negative) { most_negative = multiplier; negative = a; }
                for (int v = 0; v < variables; v++) block->delta[v] += multiplier*
                    block->gradient[(size_t)a*CS_REBUILD_BLOCK_VARS+(size_t)v]/block->norm[a];
            }
            if (negative < 0) break;
            block->active[negative] = 0;
        }
    }
    double length2 = 0;
    for (int v = 0; v < variables; v++) length2 += block->delta[v]*block->delta[v];
    if (!(length2 > 0)) for (int v = 0; v < nodes; v++) for (int d = 0; d < 3; d++)
        block->delta[3*v+d] = (double)target[3*(size_t)block->nodes[v]+(size_t)d]-block->saved[3*v+d];
    for (int v = 0; v < nodes; v++) {
        double step2 = 0;
        for (int d = 0; d < 3; d++) step2 += block->delta[3*v+d]*block->delta[3*v+d];
        if (block->trust2[v] > 0 && step2 > block->trust2[v]) alpha = fmin(alpha,sqrt(block->trust2[v]/step2));
    }
    for (int r = 0; r < faces; r++) if (block->area[r] > 0) {
        size_t f = (size_t)block->faces[r]; const int32_t *face = proposal->face_nodes+3*f;
        double delta[3][3] = {{0}}, linear = 0, quadratic = 0, e[3], t[3];
        for (int corner = 0; corner < 3; corner++) {
            int at = cs_rebuild_block_index(block->nodes,nodes,face[corner]);
            if (at >= 0) for (int d = 0; d < 3; d++) delta[corner][d] = block->delta[3*at+d];
        }
        for (int v = 0; v < variables; v++) linear += block->gradient[(size_t)r*CS_REBUILD_BLOCK_VARS+(size_t)v]*block->delta[v];
        for (int d = 0; d < 3; d++) { e[d] = delta[1][d]-delta[0][d]; t[d] = delta[2][d]-delta[0][d]; }
        for (int d = 0; d < 3; d++) quadratic += normal[3*f+(size_t)d]*(e[(d+1)%3]*t[(d+2)%3]-e[(d+2)%3]*t[(d+1)%3]);
        double zero = cs_rebuild_area_first_zero(block->area[r],linear,quadratic);
        if (zero <= alpha) alpha = .9*zero;
    }
    for (int attempt = 0; attempt < CLOSEUP_REBUILD_PROJECTION_SWEEPS && alpha > 0 && isfinite(alpha); attempt++, alpha *= .5) {
        size_t after_bad = 0, moved = 0; double after_energy = 0; int preserve = 1;
        for (int v = 0; v < nodes; v++) {
            size_t at = (size_t)block->nodes[v];
            for (int d = 0; d < 3; d++) {
                g->pos[3*at+(size_t)d] = (float)((double)block->saved[3*v+d]+alpha*block->delta[3*v+d]);
                if (!isfinite(g->pos[3*at+(size_t)d])) preserve = 0;
            }
            moved += memcmp(g->pos+3*at,block->saved+3*v,3*sizeof(float)) != 0;
        }
        for (int r = 0; r < faces; r++) {
            size_t f = (size_t)block->faces[r]; double a = cs_rebuild_motion_area(g,proposal->face_nodes+3*f,normal+3*f);
            if (!isfinite(a) || (block->area[r] > 0 && !(a > 0))) preserve = 0;
            if (!(a > 0)) { after_bad++; after_energy += a*a; }
        }
        if (preserve && moved && after_bad <= bad &&
            (after_energy < energy || (after_energy == energy && after_bad < bad))) {
            st->coupled_blocks++; st->coupled_vertices += moved; rc = 1; goto done;
        }
        for (int v = 0; v < nodes; v++) memcpy(g->pos+3*(size_t)block->nodes[v],block->saved+3*v,3*sizeof(float));
    }
done:
    Sparse_factor_free(&factor); return rc;
}
