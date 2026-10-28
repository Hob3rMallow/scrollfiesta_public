/* The field projection has one immutable matrix and eight right-hand sides.
 * Its factor belongs to this explicit scope, never to SLIM's changing matrix
 * or to scratch restored by the caller between linear solves. */
typedef struct IcProjectionLinear {
    Arena_T storage;
    SparseFactor_T factor;
    size_t nv, ne;
    double *diagonal, *blocks;
    int32_t *edge_a, *edge_b;
    int active;
    size_t phases, factorizations, solves;
} IcProjectionLinear;

static int ic_projection_begin(void *context)
{
    IcProjectionLinear *c = context;
    if (!c || c->active || c->factor || c->storage) return -1;
    /* A raising allocation precedes scope acquisition, so it owns nothing. */
    Arena_T storage = Arena_new_sized(4096);
    c->storage = storage;
    c->active = 1; c->phases++;
    return 0;
}

static void ic_projection_end(void *context)
{
    IcProjectionLinear *c = context;
    if (!c) return;
    Sparse_factor_free(&c->factor);
    if (c->storage) Arena_dispose(&c->storage);
    c->diagonal = c->blocks = NULL; c->edge_a = c->edge_b = NULL;
    c->nv = c->ne = 0; c->active = 0;
}

static int ic_projection_same(const IcProjectionLinear *c, size_t nv, size_t ne,
                               const double *diagonal, const int32_t *edge_a,
                               const int32_t *edge_b, const double *blocks)
{
    /* Pointer or hash equality cannot prove an unchanged numerical matrix. */
    return nv == c->nv && ne == c->ne &&
        !memcmp(diagonal,c->diagonal,3*nv*sizeof *diagonal) &&
        (!ne || (!memcmp(edge_a,c->edge_a,ne*sizeof *edge_a) &&
                 !memcmp(edge_b,c->edge_b,ne*sizeof *edge_b) &&
                 !memcmp(blocks,c->blocks,3*ne*sizeof *blocks)));
}

static int ic_projection_linear(Arena_T arena, void *context, size_t nv, size_t ne,
                                 const double *diagonal, const int32_t *edge_a,
                                 const int32_t *edge_b, const double *blocks,
                                 const double *rhs, double *result)
{
    IcProjectionLinear *c = context;
    if (!c || !c->active)
        return ic_slim_linear(arena,NULL,nv,ne,diagonal,edge_a,edge_b,blocks,rhs,result);
    if (!nv || nv > (size_t)INT_MAX/3 || ne > ((size_t)INT_MAX-3*nv)/4) return -1;
    /* The tiny-system fast path has different arithmetic from TAUCS. Keep it
     * byte-exact and avoid paying sparse-factor overhead for these systems. */
    if (2*nv <= (size_t)Sparse_spd_dense_limit()) {
        if (c->factor) return -1;
        c->factorizations++; c->solves++;
        return ic_slim_linear(arena,NULL,nv,ne,diagonal,edge_a,edge_b,blocks,rhs,result);
    }
    if (c->factor) {
        if (!ic_projection_same(c,nv,ne,diagonal,edge_a,edge_b,blocks)) return -1;
    } else {
        if (!c->storage) return -1;
        c->nv = nv; c->ne = ne;
        c->diagonal = ARENA_ALLOC(c->storage,3*nv*sizeof *diagonal);
        c->edge_a = ARENA_ALLOC(c->storage,(ne ? ne : 1)*sizeof *edge_a);
        c->edge_b = ARENA_ALLOC(c->storage,(ne ? ne : 1)*sizeof *edge_b);
        c->blocks = ARENA_ALLOC(c->storage,3*(ne ? ne : 1)*sizeof *blocks);
        memcpy(c->diagonal,diagonal,3*nv*sizeof *diagonal);
        if (ne) {
            memcpy(c->edge_a,edge_a,ne*sizeof *edge_a);
            memcpy(c->edge_b,edge_b,ne*sizeof *edge_b);
            memcpy(c->blocks,blocks,3*ne*sizeof *blocks);
        }
        Arena_Mark mark = Arena_save(arena);
        volatile int rc = -1;
        TRY {
            size_t count = 3*nv+4*ne, at = 0;
            int *row = ARENA_ALLOC(arena,count*sizeof *row);
            int *col = ARENA_ALLOC(arena,count*sizeof *col);
            double *value = ARENA_ALLOC(arena,count*sizeof *value);
            /* Keep the one-shot solver's exact triplet order and zero entries. */
            for (size_t v = 0; v < nv; v++) {
                int x = (int)(2*v), y = x+1;
                row[at] = col[at] = x; value[at++] = diagonal[3*v];
                row[at] = y; col[at] = x; value[at++] = diagonal[3*v+1];
                row[at] = col[at] = y; value[at++] = diagonal[3*v+2];
            }
            int valid = 1;
            for (size_t e = 0; e < ne; e++) {
                if (edge_a[e] < 0 || edge_b[e] < 0 || (size_t)edge_a[e] >= nv ||
                    (size_t)edge_b[e] >= nv || edge_a[e] == edge_b[e]) { valid = 0; break; }
                for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
                    int i = 2*edge_a[e]+a, j = 2*edge_b[e]+b;
                    row[at] = i > j ? i : j; col[at] = i > j ? j : i;
                    value[at++] = blocks[3*e+(size_t)(a+b)];
                }
            }
            if (valid) {
                c->factorizations++;
                rc = Sparse_factor_spd((int)(2*nv),(int)count,row,col,value,&c->factor);
            }
        } FINALLY { Arena_restore(arena,mark); } END_TRY;
        if (rc) return rc;
    }
    c->solves++;
    return Sparse_factor_solve(c->factor,rhs,result);
}

static void ic_projection_options(CloseupAlignOpts *options, IcProjectionLinear *cache)
{
    options->slim.linear_solve = ic_projection_linear;
    options->slim.linear_context = cache;
    options->projection_begin = ic_projection_begin;
    options->projection_end = ic_projection_end;
}
