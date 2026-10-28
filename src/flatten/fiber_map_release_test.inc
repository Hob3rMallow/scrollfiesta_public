/* Refining a shared map must add freedom, not samples, weights or canvas.
 * Test a non-affine coarse field throughout its domain, including its border,
 * against its nested fine representation and the original target objective. */
static void fm_test_value(const FiberMapResult *map, double x, double y, double out[2])
{
    int ix=(int)floor(x),iy=(int)floor(y);
    if(ix==map->cols-1)ix--;if(iy==map->rows-1)iy--;
    FmSample sample={(size_t)iy*(size_t)map->cols+(size_t)ix,x-ix,y-iy};
    fm_value(map,&sample,out);
}

static int fm_refine_selftest(Arena_T arena)
{
    enum {ROWS=3,COLS=4,N=ROWS*COLS};
    float u[N],v[N],old_u[N],old_v[N];int failures=0,pad=FIBER_MAP_PAD_CELLS;
    for(int q=0;q<N;q++){u[q]=(float)(q+3);v[q]=(float)(2*q-5);}
    u[1]=v[1]=NAN;v[6]=NAN;
    memcpy(old_u,u,sizeof u);memcpy(old_v,v,sizeof v);
    FiberMapGrid grid={ROWS,COLS,32,-48,8,u,v,{1,.125,7,0,1,-3},{2,.5}};
    for(int factor=1;factor<=4;factor*=2){
        Arena_Mark mark=Arena_save(arena);FiberMapGrid fine={0};
        if(fm_refine_grid(arena,&grid,factor,&fine)){failures++;Arena_restore(arena,mark);continue;}
        FiberMapResult a={ROWS+2*pad,COLS+2*pad,grid.row0-pad*grid.cell,grid.col0-pad*grid.cell,grid.cell,NULL};
        FiberMapResult b={fine.rows+2*pad,fine.cols+2*pad,fine.row0-pad*fine.cell,fine.col0-pad*fine.cell,fine.cell,NULL};
        if(b.rows!=factor*(a.rows-1)+1 || b.cols!=factor*(a.cols-1)+1 ||
           b.row0!=a.row0 || b.col0!=a.col0 ||
           b.row0+(b.rows-1)*b.cell!=a.row0+(a.rows-1)*a.cell ||
           b.col0+(b.cols-1)*b.cell!=a.col0+(a.cols-1)*a.cell ||
           memcmp(grid.affine,fine.affine,sizeof grid.affine) ||
           memcmp(grid.source_scale,fine.source_scale,sizeof grid.source_scale))failures++;
        size_t observed=0;
        for(int y=0;y<fine.rows;y++)for(int x=0;x<fine.cols;x++){
            size_t q=(size_t)y*(size_t)fine.cols+(size_t)x;
            if(isfinite(fine.u[q]) && isfinite(fine.v[q]))observed++;
            int dx=x-pad*(factor-1),dy=y-pad*(factor-1);
            int original=dx>=0 && dy>=0 && dx%factor==0 && dy%factor==0 && dx/factor<COLS && dy/factor<ROWS;
            if(!original){if(!isnan(fine.u[q]) || !isnan(fine.v[q]))failures++;}
            else{
                size_t p=(size_t)(dy/factor)*COLS+(size_t)(dx/factor);
                if(memcmp(fine.u+q,u+p,sizeof(float)) || memcmp(fine.v+q,v+p,sizeof(float)))failures++;
                if(fine.col0+x*fine.cell!=grid.col0+(dx/factor)*grid.cell ||
                   fine.row0+y*fine.cell!=grid.row0+(dy/factor)*grid.cell)failures++;
            }
        }
        if(observed!=N-2)failures++;
        a.uv=ARENA_ALLOC(arena,2*(size_t)a.rows*(size_t)a.cols*sizeof *a.uv);
        b.uv=ARENA_ALLOC(arena,2*(size_t)b.rows*(size_t)b.cols*sizeof *b.uv);
        for(int y=0;y<a.rows;y++)for(int x=0;x<a.cols;x++){
            size_t q=(size_t)y*(size_t)a.cols+(size_t)x;
            a.uv[2*q]=x+.3*sin(x+2*y);a.uv[2*q+1]=y+.2*cos(3*x-y);
        }
        for(int y=0;y<b.rows;y++)for(int x=0;x<b.cols;x++)
            fm_test_value(&a,(double)x/factor,(double)y/factor,b.uv+2*((size_t)y*(size_t)b.cols+(size_t)x));
        for(int y=0;y<=29;y++)for(int x=0;x<=31;x++){
            double px=(a.cols-1)*(double)x/31,py=(a.rows-1)*(double)y/29,p[2],q[2];
            fm_test_value(&a,px,py,p);fm_test_value(&b,factor*px,factor*py,q);
            if(hypot(p[0]-q[0],p[1]-q[1])>1e-12)failures++;
        }
        if(fabs(fm_target_error(&grid,&a)-fm_target_error(&fine,&b))>1e-10)failures++;
        Arena_restore(arena,mark);
    }
    FiberMapGrid invalid=grid,fine;
    invalid.rows=INT_MAX;
    if(!fm_refine_grid(arena,&invalid,2,&fine) || !fm_refine_grid(arena,&grid,0,&fine) ||
       !fm_refine_grid(arena,&grid,INT_MAX,&fine))failures++;
    invalid=grid;invalid.u=NULL;
    if(!fm_refine_grid(arena,&invalid,2,&fine) || memcmp(u,old_u,sizeof u) || memcmp(v,old_v,sizeof v))failures++;
    fprintf(stderr,"[fiber_map] nested map basis and observation preservation: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}

/* Independent regions must each release their nearest retained constraint;
 * fixing a shallow defect must not prevent a distant free region growing. */
static int fm_release_selftest(Arena_T arena)
{
    enum { ROWS=5, COLS=41, N=ROWS*COLS };
    uint8_t fixed[N]={0},seed[N]={0};FiberMapResult map={0};
    map.rows=ROWS;map.cols=COLS;int failures=0;
    for(int y=0;y<ROWS;y++)for(int x=0;x<COLS;x++)fixed[y*COLS+x]=(uint8_t)(x<=4 || x>=35);
    seed[2*COLS+4]=seed[2*COLS+25]=1;
    size_t released=fm_release(arena,&map,fixed,seed);
    if(!released || fixed[2*COLS+4] || fixed[2*COLS+35] || fixed[2*COLS+36] ||
       !fixed[2*COLS] || !fixed[2*COLS+40]){
        fprintf(stderr,"[fiber_map] release FAIL: shallow defect starves distant free region\n");failures++;
    }
    /* One central failure has two equally close boundaries. Both must grow. */
    memset(seed,0,sizeof seed);
    for(int y=0;y<ROWS;y++)for(int x=0;x<COLS;x++)fixed[y*COLS+x]=(uint8_t)(x<=5 || x>=35);
    seed[2*COLS+20]=1;
    if(!fm_release(arena,&map,fixed,seed) || fixed[2*COLS+5] || fixed[2*COLS+35] ||
       !fixed[2*COLS] || !fixed[2*COLS+40])failures++;
    uint8_t before[N];memcpy(before,fixed,sizeof fixed);memset(seed,0,sizeof seed);
    if(fm_release(arena,&map,fixed,seed) || memcmp(before,fixed,sizeof fixed))failures++;
    memset(fixed,0,sizeof fixed);seed[2*COLS+20]=1;
    if(fm_release(arena,&map,fixed,seed))failures++;
    /* Source-derived baseline pins survive releases and do not hide the
     * nearest releasable target on the far side of a pinned region. */
    memset(fixed,0,sizeof fixed);memset(seed,0,sizeof seed);
    fixed[2*COLS+20]=2;fixed[2*COLS+35]=1;seed[2*COLS+20]=1;
    if(fm_release(arena,&map,fixed,seed)!=1 || fixed[2*COLS+20]!=2 || fixed[2*COLS+35])failures++;
    fprintf(stderr,"[fiber_map] independent-region release: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}
