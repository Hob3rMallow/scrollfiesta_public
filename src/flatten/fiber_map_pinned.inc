/* A conservative recovery after the full constrained fit exhausts its budget.
 * Every knot supporting a source triangle outside the physical band is held
 * at the certified baseline. Those triangles therefore keep their original
 * metric exactly; surrounding targets may still supply useful nonrigid fit.
 * The same whole-map and original-face checks govern acceptance. Pins are
 * selected from source measurements, never from chart IDs or curated regions.
 * This is not a repair of the source distortion or a fibre-angle certificate. */
static int fm_pinned_fit(Arena_T arena,const FiberMapGrid *grid,FiberMapResult *map,
                         const double *base,const double *safe,const FmSample *sample,
                         const FmMetric *metric,const uint8_t *reversed,const int32_t *faces,
                         size_t nv,size_t nf,double *output,FiberMapStats *stats)
{
    Arena_Mark mark=Arena_save(arena);int rc=-1,pad=FIBER_MAP_PAD_CELLS;
    size_t nn=(size_t)map->rows*(size_t)map->cols,targets=0,pins=0,released=0;
    double translation[2]={0,0};size_t observed_pins=0;
    uint8_t *fixed=ARENA_CALLOC(arena,nn,1),*seed=ARENA_CALLOC(arena,nn,1);
    memcpy(map->uv,base,2*nn*sizeof *base);
    for(int y=0;y<grid->rows;y++)for(int x=0;x<grid->cols;x++){
        size_t p=(size_t)y*(size_t)grid->cols+(size_t)x;
        if(!isfinite(grid->u[p]) || !isfinite(grid->v[p]))continue;
        size_t q=(size_t)(y+pad)*(size_t)map->cols+(size_t)(x+pad);
        fixed[q]=1;targets++;map->uv[2*q]=grid->u[p];map->uv[2*q+1]=grid->v[p];
    }
    for(size_t f=0;f<nf;f++){
        if(metric[f].lower>=FIBER_MAP_METRIC_MIN-ASM_METRIC_ROUNDOFF &&
           metric[f].upper<=FIBER_MAP_METRIC_MAX+ASM_METRIC_ROUNDOFF)continue;
        for(int k=0;k<3;k++){
            size_t a=sample[(size_t)faces[3*f+(size_t)k]].node;
            size_t nodes[4]={a,a+1,a+(size_t)map->cols,a+(size_t)map->cols+1};
            for(int j=0;j<4;j++){
                size_t q=nodes[j];
                if(fixed[q]==2)continue;
                if(fixed[q]==1){
                    translation[0]+=base[2*q]-map->uv[2*q];
                    translation[1]+=base[2*q+1]-map->uv[2*q+1];observed_pins++;
                }
                released+=(size_t)(fixed[q]==1);fixed[q]=2;pins++;
                map->uv[2*q]=base[2*q];map->uv[2*q+1]=base[2*q+1];
            }
        }
    }
    if(!pins || released==targets)goto done;
    /* Target positions have an arbitrary translation, which cannot align
     * fibres. Choose that gauge from the observed baseline pins before
     * extending; otherwise an irrelevant offset creates strain at every
     * pinned boundary and can force all useful targets to be released. */
    if(observed_pins){
        translation[0]/=(double)observed_pins;translation[1]/=(double)observed_pins;
        for(size_t q=0;q<nn;q++)if(fixed[q]==1){
            map->uv[2*q]+=translation[0];map->uv[2*q+1]+=translation[1];
        }
    }
    fprintf(stderr,"[fiber_map pinned] %zu baseline knots, %zu of %zu measured targets remain; original physical bounds unchanged\n",
            pins,targets-released,targets);fflush(stderr);
    for(int pass=0;pass<FIBER_MAP_REPAIR_PASSES;pass++){
        if(fm_extend(arena,map,base,fixed))goto done;
        for(size_t v=0;v<nv;v++)fm_value(map,sample+v,output+2*v);
        memset(seed,0,nn);size_t cells=fm_cells(map,base,seed),bad=0,orientation=0;
        for(size_t f=0;f<nf;f++){
            double lo,hi;fm_metric(metric+f,faces+3*f,output,1,1,&lo,&hi);
            int physical=!fm_metric_preserved(metric+f,lo,hi);
            int32_t indices[3]={faces[3*f],faces[3*f+(reversed[f]?2:1)],faces[3*f+(reversed[f]?1:2)]};
            int folded=!UvGuard_interval(safe,output,indices,1);
            bad+=(size_t)physical;orientation+=(size_t)folded;
            if(physical || folded)for(int k=0;k<3;k++)
                fm_seed_cell(seed,sample[(size_t)faces[3*f+(size_t)k]].node,(size_t)map->cols);
        }
        fprintf(stderr,"[fiber_map pinned] pass %d: %zu metric failures, %zu triangle paths, %zu cells; %zu targets remain\n",
                pass+1,bad,orientation,cells,targets-released);fflush(stderr);
        if(!bad && !orientation && !cells){
            if(fm_harmonic_useful(grid,map,base,stats)){
                stats->baseline_pinned_nodes=pins;stats->relaxed_nodes=released;
                stats->metric_faces_after=0;stats->passes+=pass+1;rc=0;
                double gain=stats->metric_fit_before>0?1-stats->metric_fit_after/stats->metric_fit_before:1;
                fprintf(stderr,"[fiber_map pinned] certified useful fit: gain %.9g; %zu source-derived pins\n",gain,pins);fflush(stderr);
            }
            break;
        }
        size_t count=fm_release(arena,map,fixed,seed);released+=count;
        if(!count || released==targets)break;
    }
done:
    Arena_restore(arena,mark);return rc;
}

static int fm_pinned_selftest(Arena_T arena)
{
    enum {N=17,NV=N*N,NF=2*(N-1)*(N-1)};
    Arena_Mark mark=Arena_save(arena);int failures=0,pad=FIBER_MAP_PAD_CELLS;
    float u[NV],v[NV],xyz[3*NV];double canvas[2*NV],output[2*NV];int32_t faces[3*NF];
    FmSample sample[NV];FmMetric metric[NF];uint8_t reversed[NF]={0};size_t count=0;
    FiberMapGrid grid={N,N,0,0,1,u,v,{1,0,0,0,1,0},{1,1}};
    FiberMapResult map={N+2*pad,N+2*pad,-pad,-pad,1,NULL};FiberMapStats stats={0};
    size_t nn=(size_t)map.rows*(size_t)map.cols;
    double *base=ARENA_ALLOC(arena,2*nn*sizeof *base);map.uv=ARENA_ALLOC(arena,2*nn*sizeof *map.uv);
    for(int y=0;y<map.rows;y++)for(int x=0;x<map.cols;x++){
        size_t q=(size_t)y*(size_t)map.cols+(size_t)x;base[2*q]=x-pad;base[2*q+1]=y-pad;
    }
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        size_t q=(size_t)(y*N+x),node=(size_t)(y+pad)*(size_t)map.cols+(size_t)(x+pad);
        canvas[2*q]=x;canvas[2*q+1]=y;u[q]=(float)x;v[q]=(float)(y+.05*fmax(0,x-9));
        xyz[3*q]=(float)(x<=4?2*x:x+4);xyz[3*q+1]=(float)y;xyz[3*q+2]=0;
        sample[q]=(FmSample){node,0,0};
        if(x+1<N && y+1<N){
            int32_t a=(int32_t)q,b=a+1,c=a+N,d=c+1;
            faces[3*count]=a;faces[3*count+1]=b;faces[3*count+2]=d;count++;
            faces[3*count]=a;faces[3*count+1]=d;faces[3*count+2]=c;count++;
        }
    }
    for(size_t f=0;f<NF;f++){
        double lo,hi;if(fm_reference(xyz,faces+3*f,metric+f)){failures++;continue;}
        fm_metric(metric+f,faces+3*f,canvas,1,1,&lo,&hi);
        metric[f].lower=fmin(.75,lo);metric[f].upper=fmax(1.25,hi);
    }
    int rc=fm_pinned_fit(arena,&grid,&map,base,canvas,sample,metric,reversed,faces,NV,NF,output,&stats);
    if(count!=NF || rc ||
       !stats.baseline_pinned_nodes || stats.metric_faces_after || stats.metric_fit_after>.01*stats.metric_fit_before)failures++;
    for(int y=0;!rc && y<N;y++)for(int x=0;x<=4;x++){
        size_t q=(size_t)(y*N+x);
        if(hypot(output[2*q]-canvas[2*q],output[2*q+1]-canvas[2*q+1])>1e-12)failures++;
    }
    if(!rc){
        double *expected=ARENA_ALLOC(arena,2*NV*sizeof *expected);memcpy(expected,output,2*NV*sizeof *expected);
        /* Integer offsets are exactly represented by the float input grid.
         * The retained map must be independent of target-coordinate gauge. */
        for(size_t q=0;q<NV;q++){u[q]+=32;v[q]-=16;}
        rc=fm_pinned_fit(arena,&grid,&map,base,canvas,sample,metric,reversed,faces,NV,NF,output,&stats);
        if(rc)failures++;
        for(size_t q=0;!rc && q<2*NV;q++)if(fabs(expected[q]-output[q])>2e-6){failures++;break;}
    }
    fprintf(stderr,"[fiber_map] source-derived baseline pins: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    Arena_restore(arena,mark);return failures;
}
