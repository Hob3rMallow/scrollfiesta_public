/* Caller-owned, certified warm starts. A displacement is reused only for a
 * byte-identical H and b. Changed quadratics start at the current map (zero
 * displacement) and reuse only unique physical-row multiplier hints. */
struct SmLocalWarm {
    Arena_T arena;const char *capture_prefix;int valid,captured;size_t n,nt,nr;int *row,*col;
    double *a,*b,*x,*lambda;SmBound *rows;
    size_t calls,same_quadratic,changed_quadratic,matched_rows,iterations,solves,fallbacks,resets;
    double seconds,max_gap,max_dual,max_violation;
};
typedef struct SmIpmKey {
    const SmBound *key;size_t old_index,old_count,new_count;
} SmIpmKey;

static uint64_t smiw_hash_bytes(uint64_t h,const void *data,size_t bytes)
{
    const unsigned char *p=data;
    for(size_t i=0;i<bytes;i++){h^=p[i];h*=UINT64_C(1099511628211);}return h;
}
static uint64_t smiw_hash(const SmBound *r,int exact)
{
    uint64_t h=UINT64_C(1469598103934665603);
    h=smiw_hash_bytes(h,&r->face,sizeof r->face);h=smiw_hash_bytes(h,&r->upper,sizeof r->upper);
    h=smiw_hash_bytes(h,r->index,sizeof r->index);
    if(exact){h=smiw_hash_bytes(h,r->g,sizeof r->g);h=smiw_hash_bytes(h,&r->rhs,sizeof r->rhs);}return h;
}
static int smiw_equal(const SmBound *a,const SmBound *b,int exact)
{
    return a->face==b->face&&a->upper==b->upper&&!memcmp(a->index,b->index,sizeof a->index)&&
        (!exact||(!memcmp(a->g,b->g,sizeof a->g)&&!memcmp(&a->rhs,&b->rhs,sizeof a->rhs)));
}
static SmIpmKey *smiw_slot(SmIpmKey *table,size_t capacity,const SmBound *r,int exact)
{
    size_t at=(size_t)smiw_hash(r,exact)&(capacity-1);
    while(table[at].key&&!smiw_equal(table[at].key,r,exact))at=(at+1)&(capacity-1);
    if(!table[at].key)table[at].key=r;return table+at;
}
static void smiw_store(SmLocalWarm *w,size_t n,size_t nt,const int *row,const int *col,
    const double *a,const double *b,size_t nr,const SmBound *rows,const double *x,const double *lambda)
{
    w->valid=0;Arena_free(w->arena);w->n=n;w->nt=nt;w->nr=nr;
    w->row=ARENA_ALLOC(w->arena,nt*sizeof *row);w->col=ARENA_ALLOC(w->arena,nt*sizeof *col);
    w->a=ARENA_ALLOC(w->arena,nt*sizeof *a);w->b=ARENA_ALLOC(w->arena,n*sizeof *b);w->x=ARENA_ALLOC(w->arena,n*sizeof *x);
    w->rows=ARENA_ALLOC(w->arena,nr*sizeof *rows);w->lambda=ARENA_ALLOC(w->arena,nr*sizeof *lambda);
    memcpy(w->row,row,nt*sizeof *row);memcpy(w->col,col,nt*sizeof *col);memcpy(w->a,a,nt*sizeof *a);
    memcpy(w->b,b,n*sizeof *b);memcpy(w->x,x,n*sizeof *x);memcpy(w->rows,rows,nr*sizeof *rows);memcpy(w->lambda,lambda,nr*sizeof *lambda);w->valid=1;
}

static int smiw_project(Arena_T arena,size_t n,size_t nt,const int *row,const int *col,
    const double *a,const double *b,size_t nr,const SmBound *rows,double *out,double **lambda,SmMetric *metric)
{
    SmLocalWarm *w=metric->ipm;if(!w||!w->arena||!nr||nr>SIZE_MAX/8/sizeof(SmIpmKey))return -1;
    double started=ves_clock_sec();w->calls++;SparseLocalGlobalStats stats={0};
    double *initial=ARENA_CALLOC(arena,n,sizeof *initial),*initial_lambda=ARENA_CALLOC(arena,nr,sizeof *initial_lambda);
    *lambda=ARENA_CALLOC(arena,nr,sizeof **lambda);
    int available=w->valid&&w->n==n;
    int same=available&&w->nt==nt&&!memcmp(w->row,row,nt*sizeof *row)&&!memcmp(w->col,col,nt*sizeof *col)&&
        !memcmp(w->a,a,nt*sizeof *a)&&!memcmp(w->b,b,n*sizeof *b);
    if(w->valid&&w->n==n) {
        if(w->nr>SIZE_MAX/4-nr)return -1;
        size_t capacity=1;while(capacity<4*(w->nr+nr)){if(capacity>SIZE_MAX/2)return -1;capacity*=2;}
        SmIpmKey *table=ARENA_CALLOC(arena,capacity,sizeof *table);
        for(size_t r=0;r<w->nr;r++)if(same||w->rows[r].face!=SIZE_MAX) {
            SmIpmKey *key=smiw_slot(table,capacity,w->rows+r,same);key->old_index=r;key->old_count++;
        }
        for(size_t r=0;r<nr;r++)if(same||rows[r].face!=SIZE_MAX)smiw_slot(table,capacity,rows+r,same)->new_count++;
        for(size_t r=0;r<nr;r++)if(same||rows[r].face!=SIZE_MAX) {
            SmIpmKey *key=smiw_slot(table,capacity,rows+r,same);
            if(key->old_count==1&&key->new_count==1){initial_lambda[r]=w->lambda[key->old_index];w->matched_rows++;}
        }
        if(same){memcpy(initial,w->x,n*sizeof *initial);w->same_quadratic++;}else w->changed_quadratic++;
    }
    SparseLocalRow *local=ARENA_CALLOC(arena,nr,sizeof *local);
    for(size_t r=0;r<nr;r++) {
        memcpy(local[r].index,rows[r].index,sizeof local[r].index);
        memcpy(local[r].value,rows[r].g,sizeof local[r].value);local[r].lower=rows[r].rhs;
    }
    SparseLocalGlobalStats lg={0};
    int status=Sparse_local_global(arena,n,nt,row,col,a,b,nr,local,
        available?initial:NULL,available?initial_lambda:NULL,20000,out,*lambda,&lg);
    stats.converged=lg.converged;stats.iterations=lg.iterations;stats.solves=lg.solves;
    stats.primal_residual=lg.primal_residual;stats.stationarity_residual=lg.stationarity_residual;stats.complementarity=lg.complementarity;
    w->iterations+=lg.iterations;w->solves+=lg.solves;
    fprintf(stderr,"[local global solve] call=%zu iterations=%zu factors=%zu penalty=%.9g changes=%zu status=%d\n",
        w->calls,lg.iterations,lg.factors,lg.penalty,lg.penalty_changes,status);
    if(!status) {
        smiw_store(w,n,nt,row,col,a,b,nr,rows,out,*lambda);
        w->max_gap=fmax(w->max_gap,stats.complementarity);w->max_dual=fmax(w->max_dual,stats.stationarity_residual);w->max_violation=fmax(w->max_violation,stats.primal_residual);
        metric->constrained_calls++;metric->rows+=nr;metric->max_violation=fmax(metric->max_violation,stats.primal_residual);
    } else w->valid=0;
    w->seconds+=ves_clock_sec()-started;
    fprintf(stderr,"[local global warm] call=%zu same=%d matched=%zu iterations=%zu fallbacks=%zu gap=%.9g dual=%.9g status=%d\n",
        w->calls,same,w->matched_rows,w->iterations,w->fallbacks,stats.complementarity,stats.stationarity_residual,status);
    return status;
}

static int smiw_selftest(void)
{
    Arena_T arena=Arena_new();int cases=0,failures=0;
    for(int t=0;t<8;t++) {
        Arena_Mark mark=Arena_save(arena);SmLocalWarm warm={0};warm.arena=Arena_new_sized(65536);SmMetric metric={0};metric.ipm=&warm;
        double d=pow(10,(double)(t-4));int row[3]={0,1,1},col[3]={0,0,1};
        double a[3]={2*d*d,-.5,3/(d*d)},b[2]={-.45*d,.6/d};SmBound c[3]={0};
        c[0].g[0]=d;c[0].rhs=.2;c[0].face=1;c[1].index[0]=1;c[1].g[0]=1/d;c[1].rhs=.3;c[1].face=2;
        c[2].g[0]=d;c[2].rhs=.1;c[2].face=3;
        double out[2],*lambda=NULL;
        for(int pass=0;pass<3;pass++) {
            if(pass==2)b[0]=-.4*d;
            int rc=smiw_project(arena,2,3,row,col,a,b,pass?3:2,c,out,&lambda,&metric);
            failures+=rc||fabs(out[0]*d-.2)>1e-6||fabs(out[1]/d-.3)>1e-6||warm.fallbacks!=0;cases++;
        }
        failures+=warm.same_quadratic!=1||warm.changed_quadratic!=1||warm.matched_rows!=5;cases++;
        Arena_dispose(&warm.arena);Arena_restore(arena,mark);
    }
    Arena_dispose(&arena);fprintf(stderr,"sparse trajectory warm starts: %d cases, %d failures\n",cases,failures);return failures?1:0;
}
