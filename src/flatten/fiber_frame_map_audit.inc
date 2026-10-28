/* Exact derivatives of the retained map, independent of its fitting method. */
#include <float.h>

static double ff_line_error(double across, double along)
{
    return fabs(remainder(atan2(across, along), M_PI)) * 180.0 / M_PI;
}

int FiberFrame_write_directions(Arena_T arena, const FiberFrame *ff, const char *path)
{
    if (!arena || !ff || !path || ff->gh<=0 || ff->gw<=0 || ff->cell<=0 ||
        !ff->real || !ff->weight || !ff->theta || !ff->pu || !ff->pv ||
        (size_t)ff->gh>SIZE_MAX/(size_t)ff->gw/sizeof(double)) return -1;
    Arena_Mark mark=Arena_save(arena);
    size_t nn=(size_t)ff->gh*(size_t)ff->gw,nw=0,count=0;
    double *weights=ARENA_ALLOC(arena,nn*sizeof *weights);
    for (size_t i=0;i<nn;i++) if (ff->real[i] && isfinite(ff->weight[i])) weights[nw++]=ff->weight[i];
    qsort(weights,nw,sizeof *weights,ff_cmp_double);
    double median=ff_pct_sorted(weights,nw,50);
    if (!(median>0) || !isfinite(median)) { Arena_restore(arena,mark); return -1; }
    FILE *stream=fopen(path,"w");
    if (!stream) { Arena_restore(arena,mark); return -1; }
    fputs("fiber-directions-v1,column,row,h_col,h_row,v_col,v_row,weight\n",stream);
    int bad=0;
    for (int y=0;y<ff->gh;y++) for (int x=0;x<ff->gw;x++) {
        size_t i=(size_t)y*(size_t)ff->gw+(size_t)x;
        if (!ff->real[i] || !isfinite(ff->weight[i]) ||
            !((double)ff->weight[i]>FIBER_FLATTEN_EVIDENCE_REL*median)) continue;
        double e1[3],e2[3],metric[3],angle=ff->theta[i];
        if (!isfinite(angle) || ff_physical_frame(ff->pu+3*i,ff->pv+3*i,e1,e2,metric)) { bad=1; continue; }
        double hc=cos(angle),hr=sin(angle),hx=metric[0]*hc+metric[1]*hr,hy=metric[2]*hr;
        double length=hypot(hx,hy);
        if (!(length>0) || !isfinite(length)) { bad=1; continue; }
        hx/=length;hy/=length;hc/=length;hr/=length;
        double vr=hx/metric[2],vc=(-hy-metric[1]*vr)/metric[0];
        double weight=metric[0]*metric[2]*(double)ff->cell*ff->cell*fmin(1,(double)ff->weight[i]/median);
        if (!isfinite(vc) || !isfinite(vr) || !(weight>0) || !isfinite(weight)) { bad=1; continue; }
        fprintf(stream,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                (double)ff->gx0+(double)x*ff->cell,(double)ff->gy0+(double)y*ff->cell,hc,hr,vc,vr,weight);
        count++;
    }
    bad |= ferror(stream)!=0 || count==0;
    if (fclose(stream)) bad=1;
    Arena_restore(arena,mark);
    fprintf(stderr,"[fiber_frame] wrote %zu physical fibre direction observations%s\n",count,bad?" (FAILED)":"");
    return bad?-1:0;
}

/* All cells touching a coordinate on an internal grid edge. The tolerance is
 * floating-point representation error in grid coordinates, not a mesh gate. */
static int ff_map_cells(double x, int nodes, int *lo, int *hi)
{
    if (!isfinite(x) || x < 0 || x > (double)(nodes - 1)) return -1;
    double edge = floor(x + .5);
    if (fabs(x - edge) <= 16 * DBL_EPSILON * fmax(1, fabs(x))) {
        int k = (int)edge;
        *lo = k > 0 ? k - 1 : 0;
        *hi = k < nodes - 1 ? k : nodes - 2;
    } else {
        *lo = *hi = (int)floor(x);
    }
    return 0;
}

static int ff_map_jacobian(const FiberMapResult *map, int row, int col,
                            double x, double y, double J[4])
{
    size_t a = (size_t)row * (size_t)map->cols + (size_t)col;
    size_t b = a + 1, c = a + (size_t)map->cols, d = c + 1;
    double fx = x - col, fy = y - row;
    for (size_t k = 0; k < 2; k++) {
        double aa = map->uv[2*a+k], bb = map->uv[2*b+k];
        double cc = map->uv[2*c+k], dd = map->uv[2*d+k];
        if (!isfinite(aa) || !isfinite(bb) || !isfinite(cc) || !isfinite(dd)) return -1;
        J[2*k] = ((1-fy)*(bb-aa) + fy*(dd-cc)) / map->cell;
        J[2*k+1] = ((1-fx)*(cc-aa) + fx*(dd-bb)) / map->cell;
        if (!isfinite(J[2*k]) || !isfinite(J[2*k+1])) return -1;
    }
    return 0;
}

int FiberFrame_audit_map(Arena_T arena, const FiberFrame *ff,
                         const FiberMapResult *map, FiberFrameMapStats *st)
{
    if (!arena || !ff || !map || !st || ff->gh <= 0 || ff->gw <= 0 || ff->cell <= 0 ||
        !ff->real || !ff->weight || !ff->theta || !ff->pu || !ff->pv ||
        map->rows < 2 || map->cols < 2 || !map->uv || !(map->cell > 0) ||
        !isfinite(map->cell) || !isfinite(map->row0) || !isfinite(map->col0) ||
        (size_t)ff->gh > SIZE_MAX / (size_t)ff->gw / sizeof(double) ||
        (size_t)map->rows > SIZE_MAX / (size_t)map->cols / 2 / sizeof(double)) return -1;
    memset(st, 0, sizeof *st);
    st->determinant_min = INFINITY; st->determinant_max = -INFINITY;
    Arena_Mark mark = Arena_save(arena);
    size_t nn = (size_t)ff->gh * (size_t)ff->gw, nw = 0, measured = 0;
    double *weight = ARENA_ALLOC(arena, nn * sizeof *weight);
    double *ho = ARENA_ALLOC(arena, nn * sizeof *ho), *hn = ARENA_ALLOC(arena, nn * sizeof *hn);
    double *vo = ARENA_ALLOC(arena, nn * sizeof *vo), *vn = ARENA_ALLOC(arena, nn * sizeof *vn);
    for (size_t i = 0; i < nn; i++) if (ff->real[i]) {
        st->real_nodes++;
        if (isfinite(ff->weight[i])) weight[nw++] = ff->weight[i];
    }
    qsort(weight, nw, sizeof *weight, ff_cmp_double);
    double median = ff_pct_sorted(weight, nw, 50);
    for (int row = 0; row < ff->gh; row++) for (int col = 0; col < ff->gw; col++) {
        size_t i = (size_t)row * (size_t)ff->gw + (size_t)col;
        if (!ff->real[i] || !isfinite(ff->weight[i]) ||
            !((double)ff->weight[i] > FIBER_FLATTEN_EVIDENCE_REL * median)) continue;
        st->evidence_nodes++;
        double e1[3], e2[3], metric[3], angle = ff->theta[i];
        if (!isfinite(angle) || ff_physical_frame(ff->pu+3*i, ff->pv+3*i, e1, e2, metric)) {
            st->invalid_nodes++; continue;
        }
        double hc = cos(angle), hr = sin(angle);
        double hx = metric[0]*hc + metric[1]*hr, hy = metric[2]*hr;
        double length = hypot(hx, hy);
        if (!(length > 0) || !isfinite(length)) { st->invalid_nodes++; continue; }
        hx /= length; hy /= length;
        double vr = hx / metric[2], vc = (-hy - metric[1]*vr) / metric[0];
        double x = ((double)ff->gx0 + (double)col*ff->cell - map->col0) / map->cell;
        double y = ((double)ff->gy0 + (double)row*ff->cell - map->row0) / map->cell;
        int xlo = 0, xhi = 0, ylo = 0, yhi = 0, missing = 0, invalid = 0, folded = 0;
        if (ff_map_cells(x, map->cols, &xlo, &xhi) || ff_map_cells(y, map->rows, &ylo, &yhi)) {
            st->missing_nodes++; continue;
        }
        double worst_h = 0, worst_v = 0;
        for (int r = ylo; r <= yhi; r++) for (int c = xlo; c <= xhi; c++) {
            double J[4];
            if (ff_map_jacobian(map, r, c, x, y, J)) { missing = 1; continue; }
            double hu = J[0]*hc + J[1]*hr, hv = J[2]*hc + J[3]*hr;
            double vu = J[0]*vc + J[1]*vr, vv = J[2]*vc + J[3]*vr;
            double determinant = J[0]*J[3] - J[1]*J[2];
            if (!isfinite(determinant) || !isfinite(hu) || !isfinite(hv) ||
                !isfinite(vu) || !isfinite(vv) || !(hypot(hu,hv) > 0) || !(hypot(vu,vv) > 0)) {
                invalid = 1; continue;
            }
            st->cell_sides++;
            st->determinant_min = fmin(st->determinant_min, determinant);
            st->determinant_max = fmax(st->determinant_max, determinant);
            folded |= determinant <= 0;
            worst_h = fmax(worst_h, ff_line_error(hv,hu));
            worst_v = fmax(worst_v, ff_line_error(vu,vv));
        }
        st->folded_nodes += (size_t)folded;
        if (missing) { st->missing_nodes++; continue; }
        if (invalid) { st->invalid_nodes++; continue; }
        ho[measured] = ff_line_error(hr,hc); hn[measured] = worst_h;
        vo[measured] = ff_line_error(vc,vr); vn[measured] = worst_v;
        measured++;
    }
    st->supported_nodes = measured;
    st->H_old_abs = ff_quant(ho,measured); st->H_new_abs = ff_quant(hn,measured);
    st->V_old_abs = ff_quant(vo,measured); st->V_new_abs = ff_quant(vn,measured);
    if (!st->cell_sides) st->determinant_min = st->determinant_max = NAN;
    Arena_restore(arena,mark);
    return 0;
}

static int ff_map_audit_selftest(void)
{
    enum { N=25, M=81 };
    Arena_T arena = Arena_new();
    FiberFrame ff = {0}; FiberFrameMapStats st = {0};
    FiberMapResult map = {9,9,0,0,4,NULL};
    int fails = 0;
    ff.gh=ff.gw=5; ff.cell=8;
    ff.real=ARENA_ALLOC(arena,N); ff.weight=ARENA_ALLOC(arena,N*sizeof *ff.weight);
    ff.theta=ARENA_ALLOC(arena,N*sizeof *ff.theta);
    ff.pu=ARENA_CALLOC(arena,3*N,sizeof *ff.pu); ff.pv=ARENA_CALLOC(arena,3*N,sizeof *ff.pv);
    map.uv=ARENA_ALLOC(arena,2*M*sizeof *map.uv);
    for (size_t i=0;i<N;i++) { ff.real[i]=1; ff.weight[i]=1; ff.theta[i]=0; ff.pu[3*i]=1; ff.pv[3*i+1]=1; }
    /* Coarse evidence nodes see an identity map, but the retained fine cells
     * have alternating +/-.5 vertical slopes. Centered differences miss it. */
    for (int y=0;y<9;y++) for (int x=0;x<9;x++) {
        size_t i=(size_t)y*9+(size_t)x;
        map.uv[2*i]=4*x; map.uv[2*i+1]=4*y+(x%2?2:0);
    }
    int rc=FiberFrame_audit_map(arena,&ff,&map,&st);
    double expected=atan(.5)*180/M_PI;
    fails += rc || st.supported_nodes!=N || st.cell_sides!=64 || st.folded_nodes ||
             fabs(st.H_old_abs.p99)>1e-12 || fabs(st.H_new_abs.p50-expected)>1e-10 ||
             fabs(st.H_new_abs.p99-expected)>1e-10 || fabs(st.V_new_abs.p99)>1e-12;
    fprintf(stderr,"[selftest] fiber retained map fine-grid slopes: %.9g deg, expected %.9g; supported %zu\n",
            st.H_new_abs.p50,expected,st.supported_nodes);
    /* Known physical fibre axes on a sheared, stretched canvas, with large
     * input/output translations and a 180-degree output rotation. */
    ff.gx0=10000; ff.gy0=-10000; map.col0=ff.gx0; map.row0=ff.gy0;
    double angle=.23, co=cos(angle), si=sin(angle), a=2,b=.4,d=.5;
    double hr=si/d,hc=(co-b*hr)/a;
    for (size_t i=0;i<N;i++) { ff.theta[i]=(float)atan2(hr,hc); ff.pu[3*i]=a; ff.pv[3*i]=b; ff.pv[3*i+1]=d; }
    for (int y=0;y<9;y++) for (int x=0;x<9;x++) {
        size_t i=(size_t)y*9+(size_t)x;
        double wx=4*(a*x+b*y),wy=4*d*y;
        map.uv[2*i]=1e6-(co*wx+si*wy); map.uv[2*i+1]=-1e6-(-si*wx+co*wy);
    }
    rc=FiberFrame_audit_map(arena,&ff,&map,&st);
    fails += rc || st.supported_nodes!=N || st.folded_nodes || st.H_old_abs.p50<10 ||
             st.H_new_abs.p99>1e-5 || st.V_new_abs.p99>1e-5;
    /* Export the measured arms, then independently check them in physical
     * space. In a sheared canvas the vertical arm is not (-h_row,h_col). */
    {
        char path[1024];snprintf(path,sizeof path,VES_SELFTEST_DIR "/fiber_direction_export_%d.csv",ves_getpid());
        int before=fails;fails+=FiberFrame_write_directions(arena,&ff,path)!=0;
        FILE *stream=fopen(path,"r");char line[512];size_t count=0;
        if(!stream)fails++;
        else{
            if(!fgets(line,sizeof line,stream) || strcmp(line,"fiber-directions-v1,column,row,h_col,h_row,v_col,v_row,weight\n"))fails++;
            while(fgets(line,sizeof line,stream)){
                double col=0,row=0,h[2]={0},v[2]={0},weight=0;
                if(sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&col,&row,h,h+1,v,v+1,&weight)!=7){fails++;break;}
                double hx=a*h[0]+b*h[1],hy=d*h[1],vx=a*v[0]+b*v[1],vy=d*v[1];
                if(fabs(hx-co)>1e-7 || fabs(hy-si)>1e-7 || fabs(vx+si)>1e-7 || fabs(vy-co)>1e-7 ||
                   fabs(hypot(hx,hy)-1)>1e-12 || fabs(hx*vx+hy*vy)>1e-12 ||
                   weight!=a*d*64 || col!=ff.gx0+8*(double)(count%5) || row!=ff.gy0+8*(double)(count/5))fails++;
                count++;
            }
            if(ferror(stream) || count!=N)fails++;
            if(fclose(stream))fails++;
        }
        float original=ff.theta[0];ff.theta[0]=NAN;
        if(!FiberFrame_write_directions(arena,&ff,path))fails++;
        ff.theta[0]=original;
        if(fails==before)remove(path);
        fprintf(stderr,"[selftest] fiber physical direction export: %s (%d failures)\n",fails==before?"PASS":"FAIL",fails-before);
    }
    /* Missing data must remain visible instead of looking perfectly aligned. */
    map.uv[2*40]=NAN;
    rc=FiberFrame_audit_map(arena,&ff,&map,&st);
    fails += rc || st.missing_nodes!=1 || st.supported_nodes!=N-1;
    map.col0+=1000;
    rc=FiberFrame_audit_map(arena,&ff,&map,&st);
    fails += rc || st.missing_nodes!=N || st.supported_nodes || !isnan(st.H_new_abs.p50);
    map.col0=ff.gx0;
    for (int y=0;y<9;y++) for (int x=0;x<9;x++) {
        size_t i=(size_t)y*9+(size_t)x;
        map.uv[2*i]=-4*x; map.uv[2*i+1]=4*y;
    }
    rc=FiberFrame_audit_map(arena,&ff,&map,&st);
    fails += rc || st.supported_nodes!=N || st.folded_nodes!=N || !(st.determinant_max<0);
    memset(map.uv,0,2*M*sizeof *map.uv);
    rc=FiberFrame_audit_map(arena,&ff,&map,&st);
    fails += rc || st.invalid_nodes!=N || st.supported_nodes || !isnan(st.H_new_abs.p50);
    Arena_dispose(&arena);
    fprintf(stderr,"[selftest] fiber retained map audit: %s (%d failures)\n",fails?"FAIL":"PASS",fails);
    return fails;
}
