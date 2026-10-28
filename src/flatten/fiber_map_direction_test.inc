static int fm_direction_reader_selftest(Arena_T arena)
{
    const char *header="fiber-directions-v1,column,row,h_col,h_row,v_col,v_row,weight\n";
    const char *cases[]={"10,-20,1,0,0,1,4\n", "", "0,0,1,0,0,1\n", "0,0,1,0,0,-1,1\n",
        "0,0,nan,0,0,1,1\n", "0,0,1,0,0,1,0\n", "0,0,1,0,0,1,1,trailing\n",
        "0,0,1e308,0,0,1e308,1\n"};
    char path[1024];snprintf(path,sizeof path,VES_SELFTEST_DIR "/fiber_directions_%d.csv",ves_getpid());
    int failures=0;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];i++){
        Arena_Mark mark=Arena_save(arena);FILE *stream=fopen(path,"w");
        if(!stream){failures++;break;}
        fputs(header,stream);fputs(cases[i],stream);int closed=fclose(stream);
        FiberMapDirections got={0};int rc=FiberMap_read_directions(arena,path,&got);
        if(closed || (i==0?(rc || got.count!=1 || got.sample[0].column!=10 || got.sample[0].row!=-20 || got.sample[0].weight!=4):!rc))failures++;
        Arena_restore(arena,mark);
    }
    if(!failures)remove(path);
    fprintf(stderr,"[fiber_map] physical direction reader: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}

static int fm_direction_selftest(Arena_T arena)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};
    int failures=fm_direction_reader_selftest(arena);
    /* Relabeling H,V as V,-H must commute with the local minimizer. Check
     * nonorthogonal and unequal-length inputs as well as near-axis inputs;
     * signed permutations, not independent absolute values, preserve chirality. */
    for(int trial=0;trial<12;trial++){
        double a=trial==11?atan(1):.11+.09*trial,co=cos(a),si=sin(a);
        double sx=trial==11?1:1.1,sy=trial==11?1:.8,shear=trial==11?0:.13;
        double j[4]={sx*co,shear*co-sy*si,sx*si,shear*si+sy*co},z[4],q[4],p[4];
        fm_direction_prox(j,.01,1,z);
        for(int turn=0;turn<4;turn++){
            q[0]=j[1];q[1]=-j[0];q[2]=j[3];q[3]=-j[2];memcpy(j,q,sizeof j);
            q[0]=z[1];q[1]=-z[0];q[2]=z[3];q[3]=-z[2];memcpy(z,q,sizeof z);
            fm_direction_prox(j,.01,1,p);
            for(int k=0;k<4;k++)if(fabs(p[k]-z[k])>1e-13)failures++;
        }
    }
    /* Non-identical U/V systems must have independently computed residuals. */
    {
        const int row[3]={0,1,1},col[3]={0,0,1};
        const double u[3]={2,-1,2},v[3]={4,1,3},rhs[4]={0},delta[4]={1,2,1,2};
        const double expected[4]={0,-3,-6,-7};double residual[4],comp[4];
        fm_linear_residual(2,3,row,col,u,v,rhs,delta,residual,comp);
        for(int k=0;k<4;k++)if(fabs(residual[k]-expected[k])>1e-14)failures++;
    }
    /* Independent analytic checkerboard: coarse node positions alone are an
     * identity map, but every incident fine cell has a 26.565-degree H error. */
    {
        Arena_Mark mark=Arena_save(arena);double uv[2*NV];FiberMapDirection observations[25];
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            size_t q=(size_t)(y*N+x);uv[2*q]=4*x;uv[2*q+1]=4*y+(x%2?2:0);
        }
        for(int y=0;y<5;y++)for(int x=0;x<5;x++)
            observations[y*5+x]=(FiberMapDirection){8*x,8*y,{1,0},{0,1},64};
        FiberMapDirections directions={25,observations};FiberMapResult map={N,N,0,0,4,uv};
        FmDirectionStencil *stencil=NULL;size_t count=0;
        int rc=fm_direction_stencils(arena,&directions,&map,uv,&stencil,&count);
        double expected=.5*atan(.5)*atan(.5),error=rc?INFINITY:fm_direction_error(stencil,count,uv);
        if(rc || count!=64 || fabs(error-expected)>1e-13)failures++;
        for(size_t q=0;q<NV;q++){uv[2*q]=1e6-uv[2*q];uv[2*q+1]=-1e6-uv[2*q+1];}
        if(!rc && fabs(fm_direction_error(stencil,count,uv)-expected)>1e-13)failures++;
        /* The two arms cannot independently claim the same axis. */
        for(size_t q=0;q<NV;q++){uv[2*q]=4*(double)(q%N)+4*(double)(q/N);uv[2*q+1]=1e-5*4*(double)(q/N);}
        if(!rc && fm_direction_error(stencil,count,uv)<1)failures++;
        observations[0].column=-1;
        if(!fm_direction_stencils(arena,&directions,&map,uv,&stencil,&count))failures++;
        Arena_restore(arena,mark);
    }
    double unscaled_angle=INFINITY,reference[2*NV]={0};
    for(int trial=0;trial<7;trial++){
        Arena_Mark mark=Arena_save(arena);float u[NV],v[NV],xyz[3*NV];
        double canvas[2*NV],output[2*NV];int32_t faces[3*NF];FiberMapDirection observation[NV];
        double scale=trial==3 || trial==4?2:1,col0=trial==3 || trial==4?10000:0,row0=-col0;size_t nf=0;
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            size_t q=(size_t)(y*N+x);double angle=trial==2 || trial==6?0:.2;
            if(trial==1)angle+=.05*sin(.3*x)*cos(.3*y);
            double co=cos(angle),si=sin(angle);
            canvas[2*q]=col0+scale*4*x;canvas[2*q+1]=row0+scale*4*y;
            xyz[3*q]=(float)(4*x);xyz[3*q+1]=(float)(4*y);xyz[3*q+2]=0;
            /* Deliberately unaligned position target: following it cannot
             * satisfy the new physical direction objective. */
            u[q]=(float)(4*x);v[q]=(float)(4*y);
            observation[q]=(FiberMapDirection){canvas[2*q],canvas[2*q+1],{scale*co,scale*si},{-scale*si,scale*co},16};
            /* The physical cross is identical everywhere; only its arbitrary
             * local representative jumps by a quarter/half turn. */
            if(trial>=5)for(int turn=0;turn<(x+3*y)%4;turn++){
                double h[2]={observation[q].h[0],observation[q].h[1]};
                for(int d=0;d<2;d++){observation[q].h[d]=observation[q].v[d];observation[q].v[d]=-h[d];}
            }
            if(x+1<N && y+1<N){
                int32_t a=(int32_t)q,b=a+1,c=a+N,d=c+1;
                faces[3*nf]=a;faces[3*nf+1]=b;faces[3*nf+2]=d;nf++;
                faces[3*nf]=a;faces[3*nf+1]=d;faces[3*nf+2]=c;nf++;
            }
        }
        FiberMapDirections directions={NV,observation};
        FiberMapGrid grid={N,N,row0,col0,4*scale,u,v,{1/scale,0,-col0/scale,0,1/scale,-row0/scale},{scale,scale},&directions};
        FiberMapResult map={0};FiberMapStats stats={0};int rc=0;
        if(trial==4){FiberMapGrid fine;rc=fm_refine_grid(arena,&grid,2,&fine);
            if(!rc)rc=fm_fit_grid(arena,&fine,xyz,canvas,NV,faces,NF,output,&map,&stats);
        }else rc=fm_fit_grid(arena,&grid,xyz,canvas,NV,faces,NF,output,&map,&stats);
        if(rc || !stats.direction_fit || stats.direction_samples!=NV || stats.metric_faces_after ||
           !(trial==2 || trial==6?stats.direction_angle_after<1e-20:stats.direction_angle_after<.05*stats.direction_angle_before))failures++;
        if(!rc && trial!=1 && sqrt(stats.direction_angle_after)>.003)failures++;
        if(!rc && trial==0){unscaled_angle=stats.direction_angle_after;memcpy(reference,output,sizeof reference);}
        if(!rc && trial==3 && fabs(stats.direction_angle_after-unscaled_angle)>1e-8)failures++;
        if(!rc && trial==5)for(size_t q=0;q<NV;q++)for(int d=0;d<2;d++)
            if(fabs((output[2*q+(size_t)d]-output[(size_t)d])-(reference[2*q+(size_t)d]-reference[(size_t)d]))>1e-8)failures++;
        for(size_t f=0;!rc && f<nf;f++){
            FmMetric ref={0};double lo=0,hi=0;
            if(fm_reference(xyz,faces+3*f,&ref)){failures++;break;}
            ref.lower=.75;ref.upper=1.25;fm_metric(&ref,faces+3*f,output,1,1,&lo,&hi);
            if(!fm_metric_preserved(&ref,lo,hi) || !UvGuard_positive(output,faces+3*f)){failures++;break;}
        }
        fprintf(stderr,"[fiber_map] direct direction case %d: rc %d angle RMS %.9g -> %.9g rad, %zu metric failures\n",
                trial,rc,sqrt(stats.direction_angle_before),sqrt(stats.direction_angle_after),stats.metric_faces_after);
        Arena_restore(arena,mark);
    }
    fprintf(stderr,"[fiber_map] unordered quarter-turn controls: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    fprintf(stderr,"[fiber_map] physical direction fit: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}
