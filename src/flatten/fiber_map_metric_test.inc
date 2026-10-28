/* A fold in empty canvas is invisible to every original-face constraint.
 * The shared solve must restore it without suppressing the requested warp. */
static int fm_empty_cell_selftest(Arena_T arena)
{
    enum {N=9,NV=N*N,NF=8};
    Arena_Mark mark=Arena_save(arena);int failures=0,pad=FIBER_MAP_PAD_CELLS;
    float u[NV],v[NV],xyz[3*NV];double canvas[2*NV],output[2*NV];int32_t faces[3*NF];
    FmSample sample[NV];FmMetric metric[NF];uint8_t reversed[NF]={0};size_t nf=0;
    FiberMapGrid grid={N,N,0,0,1,u,v,{1,0,0,0,1,0},{1,1}};
    FiberMapResult map={N+2*pad,N+2*pad,-pad,-pad,1,NULL};FiberMapStats stats={0};
    size_t nn=(size_t)map.rows*(size_t)map.cols;
    double *base=ARENA_ALLOC(arena,2*nn*sizeof *base);
    map.uv=ARENA_ALLOC(arena,2*nn*sizeof *map.uv);
    for(int y=0;y<map.rows;y++)for(int x=0;x<map.cols;x++){
        size_t q=(size_t)y*(size_t)map.cols+(size_t)x;
        base[2*q]=x-pad;base[2*q+1]=y-pad;
    }
    memcpy(map.uv,base,2*nn*sizeof *base);
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        size_t q=(size_t)(y*N+x),node=(size_t)(y+pad)*(size_t)map.cols+(size_t)(x+pad);
        u[q]=(float)x;v[q]=(float)y;canvas[2*q]=x;canvas[2*q+1]=y;
        xyz[3*q]=(float)x;xyz[3*q+1]=(float)y;xyz[3*q+2]=0;
        sample[q]=(FmSample){node,0,0};
        if((x==0 || x==N-2) && (y==0 || y==N-2)){
            int32_t a=(int32_t)q,b=a+1,c=a+N,d=c+1;
            faces[3*nf]=a;faces[3*nf+1]=b;faces[3*nf+2]=d;nf++;
            faces[3*nf]=a;faces[3*nf+1]=d;faces[3*nf+2]=c;nf++;
        }
    }
    u[4*N+4]+=3;map.uv[2*((size_t)(4+pad)*(size_t)map.cols+(size_t)(4+pad))]+=3;
    uint8_t *seed=ARENA_CALLOC(arena,nn,1);
    if(nf!=NF || !fm_cells(&map,base,seed))failures++;
    for(size_t f=0;f<NF;f++){
        if(fm_reference(xyz,faces+3*f,metric+f))failures++;
        metric[f].lower=FIBER_MAP_METRIC_MIN;metric[f].upper=FIBER_MAP_METRIC_MAX;
    }
    int rc=fm_metric_fit(arena,&grid,&map,base,canvas,sample,metric,reversed,faces,NV,NF,output,&stats);
    if(rc || fm_cells(&map,base,seed) || stats.metric_faces_after ||
       !(stats.metric_fit_after<.75*stats.metric_fit_before))failures++;
    /* A rigid translation must not buy apparent progress toward the warp. */
    double error=fm_target_error(&grid,&map);
    for(size_t q=0;q<nn;q++){map.uv[2*q]+=137;map.uv[2*q+1]-=91;}
    if(fabs(fm_target_error(&grid,&map)-error)>1e-9)failures++;
    /* Affine cell derivatives agree at all four corners.
     * A positive half-turn endpoint still has an inadmissible path. */
    for(size_t k=0;k<4;k++){
        FmStencil s;fm_cell_stencil(k,&map,base,&s);
        double identity[4]={1,0,0,1},half[4]={-1,0,0,-1},out[4];
        for(int j=0;j<4;j++)if(fabs(s.base[j]-identity[j])>1e-12)failures++;
        if(fm_cell_path(s.base,half,1) || fm_project_cell(half,&s,out) || !fm_cell_path(s.base,out,1)){
            fprintf(stderr,"[fiber_map] half-turn cell corner %zu projection/path failed\n",k);failures++;
        }
    }
    fprintf(stderr,"[fiber_map] empty-canvas cell recovery: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    Arena_restore(arena,mark);return failures;
}

/* Check the convex projection against its defining variational inequality,
 * using independently constructed feasible matrices. This catches a merely
 * feasible clipping rule that is not the required nearest-point projection. */
static int fm_inner_selftest(void)
{
    int failures=0;
    {
        const int row[4]={0,1,1,1},col[4]={0,0,1,1};
        const double value[4]={2,-1,.5,1.5},rhs[4]={1,0,0,1},delta[4]={.2,.1,0,0};
        const double expected[4]={.7,0,0,1};double residual[4],compensation[4];
        fm_linear_residual(2,4,row,col,value,NULL,rhs,delta,residual,compensation);
        for(int k=0;k<4;k++)if(fabs(residual[k]-expected[k])>1e-14)failures++;
    }
    for(int old=0;old<2;old++)for(int reverse=0;reverse<2;reverse++){
        double lower=old?.5:.75,upper=old?2:1.25,sign=reverse?-1:1,cs=.6,sn=.8;
        double base[4]={cs*lower,-sn*upper,sign*sn*lower,sign*cs*upper},fixed[4];
        if(fm_project_inner(base,base,lower,upper,reverse,fixed)){failures++;continue;}
        for(int k=0;k<4;k++)if(fabs(fixed[k]-base[k])>1e-12)failures++;
        for(int trial=0;trial<96;trial++){
            double input[4],out[4],scaled[4],large_base[4],large[4];
            for(int k=0;k<4;k++){
                input[k]=trial==0?0:3*sin((trial+1)*(k+1)*1.713);
                scaled[k]=3*input[k];large_base[k]=3*base[k];
            }
            if(fm_project_inner(input,base,lower,upper,reverse,out) ||
               fm_project_inner(scaled,large_base,3*lower,3*upper,reverse,large)){failures++;continue;}
            for(int k=0;k<4;k++)if(fabs(large[k]-3*out[k])>1e-10)failures++;
            double a=cs*out[0]+sn*sign*out[2],b=cs*out[1]+sn*sign*out[3];
            double c=-sn*out[0]+cs*sign*out[2],d=-sn*out[1]+cs*sign*out[3];
            double eig=.5*(a+d-hypot(a-d,b+c));
            if(eig<lower-1e-10)failures++;
            for(int step=0;step<=4;step++){
                double j[4],alpha=.25*step;
                for(int k=0;k<4;k++)j[k]=base[k]+alpha*(out[k]-base[k]);
                double hi=.5*(hypot(j[0]+j[3],j[2]-j[1])+hypot(j[0]-j[3],j[2]+j[1]));
                double det=j[0]*j[3]-j[1]*j[2],lo=fabs(det)/hi;
                if(lo<lower-1e-10 || hi>upper+1e-10 || sign*det<=0)failures++;
            }
            for(int witness=0;witness<48;witness++){
                double t=.5*(upper-lower)*(witness%7)/6;
                double p=lower+t+(upper-lower-2*t)*((witness/7)%5)/4;
                double q=sqrt(fmax(0,(upper-t)*(upper-t)-p*p))*sin(witness*2.13);
                double r=t*cos(witness*1.71),s=t*sin(witness*1.71);
                double feasible[4]={cs*(p+r)-sn*(s+q),cs*(s-q)-sn*(p-r),
                    sign*(sn*(p+r)+cs*(s+q)),sign*(sn*(s-q)+cs*(p-r))};
                double dot=0;
                for(int k=0;k<4;k++)dot+=(out[k]-input[k])*(feasible[k]-out[k]);
                if(dot< -1e-9)failures++;
            }
        }
    }
    fprintf(stderr,"[fiber_map] convex polar projection: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}

static int fm_harmonic_quality_selftest(void)
{
    enum {N=2};
    int pad=FIBER_MAP_PAD_CELLS,side=N+2*FIBER_MAP_PAD_CELLS,failures=0;
    float u[4]={0,2,0,1},v[4]={0,0,1,1};
    double base[(N+2*FIBER_MAP_PAD_CELLS)*(N+2*FIBER_MAP_PAD_CELLS)*2];
    double mapped[sizeof base/sizeof *base];
    FiberMapGrid grid={N,N,0,0,1,u,v,{1,0,0,0,1,0},{1,1}};
    FiberMapResult map={side,side,-pad,-pad,1,mapped};FiberMapStats stats={0};
    for(int y=0;y<side;y++)for(int x=0;x<side;x++){
        size_t q=(size_t)(y*side+x);base[2*q]=x-pad;base[2*q+1]=y-pad;
        mapped[2*q]=base[2*q]+17;mapped[2*q+1]=base[2*q+1]-11;
    }
    failures+=fm_harmonic_useful(&grid,&map,base,&stats)!=0;
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        size_t q=(size_t)((y+pad)*side+x+pad),p=(size_t)(y*N+x);
        mapped[2*q]=u[p]+17;mapped[2*q+1]=v[p]-11;
    }
    failures+=fm_harmonic_useful(&grid,&map,base,&stats)==0;
    u[1]=1;
    for(size_t i=0;i<sizeof base/sizeof *base;i++)mapped[i]=base[i];
    failures+=fm_harmonic_useful(&grid,&map,base,&stats)==0;
    fprintf(stderr,"[fiber_map] harmonic useful-fit controls: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}

static int fm_polish_center_selftest(void)
{
    const int row[4]={0,1,1,1},col[4]={0,0,1,1};
    const double a[4]={3,-1,2,2},b[4]={2,.5,1,2},x[4]={2,3,5,7};
    const double expected[4]={3,10,13.5,23.5};
    double rhs[4]={0};FmStencil stencil[2]={{0}};int failures=0;
    for(size_t i=0;i<2;i++)for(size_t k=0;k<4;k++)stencil[i].dual[k]=(double)(i+k+1);
    fm_polish_center(2,4,row,col,a,b,x,rhs,stencil,2);
    for(size_t i=0;i<4;i++)failures+=fabs(rhs[i]-expected[i])>1e-12;
    for(size_t i=0;i<2;i++)for(size_t k=0;k<4;k++)failures+=fabs(stencil[i].dual[k])>0;
    fprintf(stderr,"[fiber_map] feasibility polish center: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}

static int fm_metric_selftest(Arena_T arena)
{
    int failures=fm_inner_selftest()+fm_empty_cell_selftest(arena);FmMetric bounds={1,0,1,.75,1.25};
    {
        /* A rigid output rotation preserves both singular values, including
         * a source face already at both of its preservation bounds. The fixed
         * polar convex subset excludes it; the actual physical set does not. */
        FmMetric reference={1,0,1,.5,2};double base[4]={.5,0,0,2};
        double co=cos(.2),si=sin(.2),rotated[4]={.5*co,-2*si,.5*si,2*co},exact[4],inner[4];
        int rc=fm_project(rotated,&reference,0,exact),ic=fm_project_inner(rotated,base,.5,2,0,inner);
        double changed=0,difference=0;
        for(int k=0;k<4;k++){changed+=fabs(inner[k]-rotated[k]);difference+=fabs(exact[k]-rotated[k]);}
        int okay=!rc && !ic && difference<1e-12 && changed>1e-3 && fm_cell_path(base,rotated,1) &&
            fabs(hypot(rotated[0],rotated[2])-.5)<1e-12 && fabs(hypot(rotated[1],rotated[3])-2)<1e-12 &&
            fabs(rotated[0]*rotated[1]+rotated[2]*rotated[3])<1e-12;
        failures+=!okay;
        fprintf(stderr,"[fiber_map] active-bound rotation: %s; physical projection change %.9g, fixed-frame change %.9g\n",
                okay?"ok":"FAIL",difference,changed);
    }
    const double cases[][4]={{.4,0,0,2},{0,-2,.4,0},{.4,0,0,-2},{0,0,0,0},{1,.2,0,1}};
    for(size_t c=0;c<sizeof cases/sizeof cases[0];c++)for(int reverse=0;reverse<2;reverse++){
        double p[4];
        if(fm_project(cases[c],&bounds,reverse,p)){failures++;continue;}
        double hi=.5*(hypot(p[0]+p[3],p[2]-p[1])+hypot(p[0]-p[3],p[2]+p[1]));
        double det=p[0]*p[3]-p[1]*p[2],lo=fabs(det)/hi;
        if(!fm_metric_preserved(&bounds,lo,hi) || ((det<0)!=reverse))failures++;
    }
    /* A projection onto the preservation set must fix every feasible input.
     * Contracting the interval asks for an unrequested strict improvement;
     * neighbouring faces sharing map coefficients may have no common interior. */
    for(int old=0;old<2;old++)for(int reverse=0;reverse<2;reverse++)for(int rotated=0;rotated<2;rotated++){
        FmMetric reference={1,0,1,old?.5:.75,old?2:1.25};
        double sign=reverse?-1:1;
        double j[4]={reference.lower,0,0,sign*reference.upper},p[4];
        if(rotated){j[0]=j[3]=0;j[1]=-reference.upper;j[2]=sign*reference.lower;}
        if(fm_project(j,&reference,reverse,p)){failures++;continue;}
        for(int k=0;k<4;k++)if(fabs(p[k]-j[k])>1e-12)failures++;
    }
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};
    float u[NV],v[NV],xyz[3*NV];double canvas[2*NV],output[2*NV];int32_t faces[3*NF];
    for(int trial=0;trial<3;trial++){
      size_t count=0;
      for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        size_t q=(size_t)(y*N+x);canvas[2*q]=4*x;canvas[2*q+1]=4*y;
        /* The second case has pre-existing 2:1 physical compression in its
         * left half. It may not buy further shrink there to satisfy fibres. */
        xyz[3*q]=(float)(trial?(x<=4?8*x:32+4*(x-4)):4*x);
        xyz[3*q+1]=(float)(4*y);xyz[3*q+2]=0;
        u[q]=(float)(.55*4*x);v[q]=(float)(4*y);
        if(x+1<N && y+1<N){
            int32_t a=(int32_t)q,b=a+1,c=a+N,d=c+1;
            faces[3*count]=a;faces[3*count+1]=b;faces[3*count+2]=d;count++;
            faces[3*count]=a;faces[3*count+1]=d;faces[3*count+2]=c;count++;
        }
      }
      FiberMapGrid grid={N,N,0,0,4,u,v,{1,0,0,0,1,0},{1,1}};
      FiberMapResult map={0};FiberMapStats stats={0};Arena_Mark mark=Arena_save(arena);
      int rc;
      if(trial==2){
          /* Always exercise the refined heterogeneous regression, including
           * when production uses the coarse basis. Its objective used to
           * stagnate before the known feasible witness's bounds were met. */
          FiberMapGrid fine;
          rc=fm_refine_grid(arena,&grid,2,&fine);
          if(!rc)rc=fm_fit_grid(arena,&fine,xyz,canvas,NV,faces,NF,output,&map,&stats);
      }else rc=FiberMap_fit(arena,&grid,xyz,canvas,NV,faces,NF,output,&map,&stats);
      /* Independent feasible witness: keep the already compressed half at
       * its original width and compress the remaining material to 0.75.
       * Compare after removing translation; the heterogeneous case cannot
       * halve that error (its witness gain is only 45.839%). */
      double differences[N],mean=0,witness_error=0;
      for(int x=0;x<N;x++){
        double known=trial?(x<=4?4*x:16+3*(x-4)):3*x;
        differences[x]=known-(double)u[x];mean+=differences[x]/N;
      }
      for(int x=0;x<N;x++)witness_error+=N*(differences[x]-mean)*(differences[x]-mean);
      if(rc || !stats.metric_iterations ||
       stats.metric_faces_after || stats.metric_fit_after>witness_error+1e-4*stats.metric_fit_before || !UvGuard_mesh(output,faces,NF)){
        fprintf(stderr,"[fiber_map] compression control %d: rc %d, metric %zu, error %.12g / %.12g\n",
                trial,rc,stats.metric_faces_after,stats.metric_fit_after,stats.metric_fit_before);failures++;
      }
      /* A pure translation cannot satisfy this improvement check. */
      if(!rc && output[2*(4*N+8)]-output[2*(4*N+4)]>=.99*16)failures++;
      for(size_t f=0;!rc && f<NF;f++){
        FmMetric reference={0};double lo,hi;
        if(fm_reference(xyz,faces+3*f,&reference)){failures++;break;}
        fm_metric(&reference,faces+3*f,canvas,1,1,&lo,&hi);
        reference.lower=fmin(.75,lo);reference.upper=fmax(1.25,hi);
        fm_metric(&reference,faces+3*f,output,1,1,&lo,&hi);
        if(!fm_metric_preserved(&reference,lo,hi)){failures++;break;}
      }
      Arena_restore(arena,mark);
    }
    fprintf(stderr,"[fiber_map] constrained original-face map fit: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}
