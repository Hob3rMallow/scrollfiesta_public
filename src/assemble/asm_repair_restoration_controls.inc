/* A metric-valid bent strip cannot fit between the fixed walls by a rigid
 * pose. Its original straight shape can. A failed current-shape forest must
 * not discard the original solver-entry alternative. */
static int ars_entry_forest_control(void)
{
    enum { NX=33,NY=9,NV=NX*NY,NF=2*(NX-1)*(NY-1) };
    Arena_T arena=Arena_new();AsmRun run={0};AsmChart charts[4]={{0}};
    float xyz[2][3*NV],uv[2][2*NV],walls[2][12],wall_uv[2][8];
    int32_t faces[3*NF],wall_faces[]={0,1,2,0,2,3};AsmCorr corr[NY]={{0}};
    for(int c=0;c<2;c++)for(int y=0;y<NY;y++)for(int x=0;x<NX;x++){
        int v=y*NX+x;double u=180.*x/(NX-1)+(c?0:-210),w=y-4.;
        xyz[c][3*v]=(float)u;xyz[c][3*v+1]=(float)w;xyz[c][3*v+2]=0;
        uv[c][2*v]=(float)u;uv[c][2*v+1]=(float)w;
        if(!c && x+1<NX && y+1<NY){int32_t tri[]={v,v+1,v+NX+1,v,v+NX+1,v+NX};memcpy(faces+6*(y*(NX-1)+x),tri,sizeof tri);}
    }
    for(int c=0;c<4;c++){
        AsmChart *ch=charts+c;ch->id=c;ch->parent=-1;ch->component=c;ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;
        if(c<2){ch->nv=NV;ch->nf=NF;ch->xyz=xyz[c];ch->uv=uv[c];ch->faces=faces;}
        else{
            int at=c-2;double low=at?-400:10,high=at?-10:400;
            double points[]={-400,low,400,low,400,high,-400,high};
            for(int v=0;v<4;v++){walls[at][3*v]=(float)points[2*v];walls[at][3*v+1]=(float)points[2*v+1];walls[at][3*v+2]=(float)c;
                wall_uv[at][2*v]=(float)points[2*v];wall_uv[at][2*v+1]=(float)points[2*v+1];}
            ch->nv=4;ch->nf=2;ch->xyz=walls[at];ch->uv=wall_uv[at];ch->faces=wall_faces;
        }
    }
    for(int i=0;i<NY;i++){
        AsmCorr *w=corr+i;w->va=i*NX+NX-1;w->vb=i*NX;w->valid=3;w->gap_a[0]=w->gap_b[0]=30;
        for(int c=0;c<2;c++){
            int v=c?w->vb:w->va;AsmTrim *trim=c?&w->trim_b:&w->trim_a;
            for(int t=0;t<NF;t++)if(faces[3*t]==v || faces[3*t+1]==v || faces[3*t+2]==v){
                const float *a=xyz[c]+3*faces[3*t],*b=xyz[c]+3*faces[3*t+1],*d=xyz[c]+3*faces[3*t+2];
                double x=b[0]-a[0],y=b[1]-a[1],u=d[0]-a[0],z=d[1]-a[1],det=x*z-y*u;
                trim->face=t;trim->valid=1;trim->inverse[0]=z/det;trim->inverse[1]=-u/det;
                trim->inverse[2]=-y/det;trim->inverse[3]=x/det;break;
            }
        }
    }
    AsmRelation relation={0};relation.a=0;relation.b=1;relation.corr_count=NY;relation.n_corr=NY;relation.continuity=ASM_CONT_SOURCE;
    run.arena=arena;run.charts=charts;run.n_charts=4;run.rels=&relation;run.n_rels=1;run.corr=corr;run.n_corr=NY;
    AsmField *f=AsmField_new(arena,&run);int fails=0;if(!f){Arena_dispose(&arena);return 1;}
    size_t bytes=2*f->nv*sizeof(double);double *initial=ARENA_ALLOC(arena,bytes),*entry=ARENA_ALLOC(arena,bytes);memcpy(initial,f->uv,bytes);
    double radius=180/3.141592653589793;
    for(int y=0;y<NY;y++)for(int x=0;x<NX;x++){
        size_t v=f->vertex_offset[1]+y*NX+x;double theta=3.141592653589793*x/(NX-1),r=radius-(y-4.);
        f->uv[2*v]=r*sin(theta);f->uv[2*v+1]=-r*cos(theta);
    }
    memcpy(entry,f->uv,bytes);uint8_t fixed_charts[]={1,0,1,1};AsmRepairStats stats={0};AprSolver s={0};s.field=f;s.stats=&stats;
    AsmRepair_defaults(&s.opts);s.opts.fixed_charts=fixed_charts;
    s.fixed=ARENA_ALLOC(arena,f->nv);s.rigid_fixed=fixed_charts;s.active_faces=ARENA_ALLOC(arena,NF*sizeof(size_t));
    uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
    for(size_t v=0;v<f->nv;v++)s.fixed[v]=(uint8_t)(f->vertex_chart[v]!=1);
    for(size_t t=0;t<f->nf;t++)if(f->face_chart[t]==1){
        s.active_faces[s.n_active++]=t;active[t]=1;double j[4],lo,hi,det;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
        if(lo<.9 || hi>1.1 || !UvGuard_positive(f->uv,f->faces+3*t))fails++;
    }
    AsmContacts_active_faces(f->contacts,active);
    if(apr_initialize_retry(&s,entry,0)!=0)fails++;
    memcpy(f->uv,entry,bytes);memset(&stats,0,sizeof stats);s.initial_uv=initial;
    if(apr_initialize_retry(&s,entry,0)!=1 || stats.forest_entry_restarts!=1 || stats.forest_initialized!=2)fails++;
    AsmSeamMeasure seam;AsmField_seam(f,0,f->uv,&seam);AsmContactStats contact={0};
    if(!seam.pass || seam.rms>1e-8 || AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&contact) || contact.pairs)fails++;
    for(size_t v=0;v<f->nv;v++)if(s.fixed[v] && memcmp(f->uv+2*v,initial+2*v,2*sizeof(double)))fails++;
    for(size_t k=0;k<s.n_active;k++){double j[4],lo,hi,det;AsmField_face(f,s.active_faces[k],f->uv,j,&lo,&hi,&det);if(lo<.9 || hi>1.1)fails++;}
    Arena_dispose(&arena);fprintf(stderr,"  native original-entry forest after obstructed metric-valid deformation: %s (%d failures)\n",fails?"FAIL":"ok",fails);return fails;
}

/* Continuous contact must catch a crossing even when both endpoint fields
 * are separated. Exercise the actual boundary hierarchy and original mesh. */
static int ars_boundary_control(void)
{
    int fails=0;
    for(int transformed=0;transformed<2;transformed++){
        Arena_T arena=Arena_new();AsmRun run={0};AsmChart charts[2]={{0}};
        float source[36],local[24];double placed[24],direction[24]={0},trial[24];
        int32_t faces[]={0,1,4,0,4,3,1,2,5,1,5,4,3,4,7,3,7,6,4,5,8,4,8,7,0,1,2};
        for(int y=0;y<3;y++)for(int x=0;x<3;x++){int v=3*y+x;local[2*v]=(float)x;local[2*v+1]=(float)y;}
        const float moving[]={-2,.5f,-1.5f,.5f,-2,1};memcpy(local+18,moving,sizeof moving);
        for(int k=0;k<24;k++)local[k]*=16;
        double angle=transformed?.43:0,cs=cos(angle),sn=sin(angle);
        for(int v=0;v<12;v++){
            source[3*v]=local[2*v];source[3*v+1]=local[2*v+1];source[3*v+2]=v<9?0.f:10.f;
            placed[2*v]=cs*local[2*v]-sn*local[2*v+1]+(transformed?10000:0);
            placed[2*v+1]=sn*local[2*v]+cs*local[2*v+1]-(transformed?3000:0);
            if(v>=9){direction[2*v]=80*cs;direction[2*v+1]=80*sn;}
        }
        for(int c=0;c<2;c++){
            charts[c].id=c;charts[c].parent=-1;charts[c].component=c;charts[c].nv=c?3:9;charts[c].nf=c?1:8;
            charts[c].xyz=source+(c?27:0);charts[c].uv=local+(c?18:0);charts[c].placed_uv=placed+(c?18:0);
            charts[c].faces=faces+(c?24:0);charts[c].placed=1;charts[c].placement_state=ASM_PLACE_ROOT;
        }
        run.arena=arena;run.charts=charts;run.n_charts=2;AsmField *f=AsmField_new(arena,&run);
        if(!f){fails++;Arena_dispose(&arena);continue;}
        uint8_t fixed[]={1,0};AprSolver s={0};s.field=f;s.opts.fixed_charts=fixed;
        ArsProblem p={0};p.arena=arena;p.s=&s;
        if(!ars_boundary_new(&p)){fails++;goto done;}
        if(!ars_boundary_energy(&p,f->uv) || fabs(p.energy)>1e-8)fails++;
        double limit=ars_boundary_limit(&p,f->uv,direction);
        if(fabs(limit-.3)>1e-9)fails++;
        for(int phase=0;phase<4;phase++){
            double alpha=phase==0?0:phase==1?1:phase==2?.5:.8*limit;
            for(int k=0;k<24;k++)trial[k]=f->uv[k]+alpha*direction[k];
            AsmContactStats contacts={0};
            int rc=AsmContacts_measure(f->contacts,trial,NULL,NULL,NULL,&contacts);
            if(rc || (phase==2?!contacts.pairs:contacts.pairs)){fails++;fprintf(stderr,"  swept boundary pose %d phase %d: rc %d, %zu pairs complete %d\n",transformed,phase,rc,contacts.pairs,contacts.complete);}
        }
        if(memcmp(f->uv,placed,sizeof placed))fails++;
done:
        if(p.boundary)free(p.boundary->edges);Arena_dispose(&arena);
    }
    fprintf(stderr,"  native swept boundary crossing with separated endpoints: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    return fails;
}

static int ars_stretch_control(void)
{
    int fails=0;const double examples[][4]={{1,0,0,1},{1.1,.07,-.03,.89},{.2,.04,-.1,1.5},{.9,-.4,.4,.9},{1.1,0,0,1.05}};
    for(size_t at=0;at<sizeof examples/sizeof *examples;at++){
        double j[4],gradient[4],hessian[16];memcpy(j,examples[at],sizeof j);ArsStretch band={0};
        if(!ars_stretch_update(&band,j,3)){fails++;continue;}
        double energy=ars_stretch_face(&band,j,1.7,gradient,hessian);if(!isfinite(energy)){fails++;continue;}
        for(int k=0;k<4;k++){
            double old=j[k],step=1e-6,gp[4],gm[4];
            j[k]=old+step;double plus=ars_stretch_face(&band,j,1.7,gp,NULL);
            j[k]=old-step;double minus=ars_stretch_face(&band,j,1.7,gm,NULL);j[k]=old;
            if(fabs((plus-minus)/(2*step)-gradient[k])>2e-6*fmax(1,fabs(gradient[k])))fails++;
            /* This tensile field has positive unprojected curvature. */
            if(at==4)for(int q=0;q<4;q++)if(fabs((gp[q]-gm[q])/(2*step)-hessian[4*q+k])>2e-6*fmax(1,fabs(hessian[4*q+k])))fails++;
        }
        for(int probe=0;probe<16;probe++){
            double square=0;for(int a=0;a<4;a++)for(int b=0;b<4;b++){
                if(fabs(hessian[4*a+b]-hessian[4*b+a])>1e-10)fails++;
                square+=((probe>>a)&1?1:-1)*hessian[4*a+b]*((probe>>b)&1?1:-1);
            }
            if(square< -1e-9 || !isfinite(square))fails++;
        }
    }
    ArsStretch bands[2]={{0}};double first[]={.6,0,0,1.4},second[]={.9,0,0,1.1};
    if(!ars_stretch_update(bands,first,1) || !ars_stretch_update(bands+1,second,1) ||
        !(bands[0].lower<bands[1].lower && bands[0].upper>bands[1].upper))fails++;
    first[0]=.8;first[3]=1.2;if(!ars_stretch_update(bands,first,1) || !bands[0].lower_qualified || !bands[0].upper_qualified)fails++;
    first[0]=.74;if(ars_stretch_update(bands,first,1) || isfinite(ars_stretch_face(bands,first,1,NULL,NULL)))fails++;
    first[0]=1;first[3]=1.26;if(ars_stretch_update(bands,first,1) || isfinite(ars_stretch_face(bands,first,1,NULL,NULL)))fails++;
    ArsPressure pressure={0};pressure.weight=1;
    for(int k=0;k<6;k++)if(ars_pressure_begin(&pressure,1))fails++;
    if(!ars_pressure_begin(&pressure,1) || pressure.weight!=4)fails++;
    for(int k=0;k<6;k++)if(ars_pressure_begin(&pressure,.5))fails++;
    if(pressure.weight!=4 || !ars_pressure_increase(&pressure,.5) || pressure.weight!=16 || ars_pressure_increase(&pressure,0))fails++;
    fprintf(stderr,"  native original-face stretch derivatives, envelopes and adaptive pressure: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    return fails;
}

static int ars_derivative_control(void)
{
    const double examples[][4]={{1,0,0,1},{1.1,.07,-.03,.89},{.2,.04,-.1,1.5},{.9,-.4,.4,.9}};
    int fails=0;
    for(size_t at=0;at<sizeof examples/sizeof *examples;at++){
        double j[4],gradient[4],hessian[16];memcpy(j,examples[at],sizeof j);
        double energy=ars_sd(j,1.7,gradient,hessian);if(!isfinite(energy)){fails++;continue;}
        for(int k=0;k<4;k++){
            double old=j[k],step=1e-6;j[k]=old+step;double plus=ars_sd(j,1.7,NULL,NULL);
            j[k]=old-step;double minus=ars_sd(j,1.7,NULL,NULL);j[k]=old;
            if(fabs((plus-minus)/(2*step)-gradient[k])>1e-6*fmax(1,fabs(gradient[k])))fails++;
        }
        for(int probe=0;probe<16;probe++){
            double square=0;
            for(int a=0;a<4;a++)for(int b=0;b<4;b++){
                if(fabs(hessian[4*a+b]-hessian[4*b+a])>1e-10)fails++;
                square+=((probe>>a)&1?1:-1)*hessian[4*a+b]*((probe>>b)&1?1:-1);
            }
            if(!(square>0) || !isfinite(square))fails++;
        }
    }
    const double inverted[]={-1,0,0,1},collapsed[]={0,0,0,1};
    if(isfinite(ars_sd(inverted,1,NULL,NULL)) || isfinite(ars_sd(collapsed,1,NULL,NULL)))fails++;
    fprintf(stderr,"  native original-area restoration derivatives and positive curvature: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    return fails;
}
