/* Packing controls protect connectivity through absent source charts, as
 * well as passing deferred joins that are not in the required graph. */
static int anp_ribbon_selftest(void)
{
    int fail=0;
    for(int test=0;test<3;test++){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.continuity_ready=1;run.layer_d_global=test==1?0:10;
        AsmChart ch[3]={{0}};float xyz[3][12],uv[3][8];int32_t tri[]={0,1,2,0,2,3};uint8_t boundary[]={1,1,1,1};
        const double radius[]={300,100,200},height[]={20,10,30};
        for(size_t c=0;c<3;c++){
            ch[c].id=(int32_t)c;ch[c].nv=4;ch[c].nf=2;ch[c].area3d=4*height[c];ch[c].xyz=xyz[c];ch[c].uv=uv[c];
            ch[c].faces=tri;ch[c].boundary=boundary;ch[c].placed=1;ch[c].placement_state=ASM_PLACE_ROOT;
            ch[c].pose_theta=.4*(double)(c+1);ch[c].pose_x=1000*(double)c;
            for(int v=0;v<4;v++){
                double x=v==1 || v==2?2:-2,z=v>=2?height[c]:0;
                xyz[c][3*v]=(float)z;xyz[c][3*v+1]=(float)sqrt(radius[c]*radius[c]-4);xyz[c][3*v+2]=(float)x;
                uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)z;
                for(int d=0;d<3;d++)ch[c].centroid[d]+=xyz[c][3*v+d]/4.0;
            }
        }
        run.charts=ch;run.n_charts=3;double origin[]={0,0,0},direction[]={1,0,0};
        AsmAxis *axis=AsmAxis_line(arena,origin,direction,0,100);
        AsmField *f=AsmField_new_local(arena,&run,NULL,NULL);AnpStats stats;
        if(!f){fail++;Arena_dispose(&arena);continue;}
        double saved[24];memcpy(saved,f->uv,sizeof saved);AsmChart old[3];memcpy(old,ch,sizeof old);
        int kept=anp_pack_layout(&run,f,8,test==2?ves_clock_sec()-1:0,&stats,NULL,1,axis,NULL);
        int okay=test?(!kept && !memcmp(saved,f->uv,sizeof saved) && !memcmp(old,ch,sizeof old)):kept;
        if(!test){
            double lo[3]={DBL_MAX,DBL_MAX,DBL_MAX},hi[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
            for(size_t c=0;c<3;c++)for(size_t v=0;v<4;v++){
                size_t at=f->vertex_offset[c]+v;lo[c]=fmin(lo[c],f->uv[2*at]);hi[c]=fmax(hi[c],f->uv[2*at]);
                for(size_t w=0;w<v;w++){
                    size_t other=f->vertex_offset[c]+w;
                    double a=hypot(saved[2*at]-saved[2*other],saved[2*at+1]-saved[2*other+1]);
                    double b=hypot(f->uv[2*at]-f->uv[2*other],f->uv[2*at+1]-f->uv[2*other+1]);
                    okay=okay && fabs(a-b)<1e-9;
                }
            }
            okay=okay && stats.regions==3 && fabs(stats.width-28)<1e-3 && fabs(stats.height-30)<1e-3 &&
                fabs(lo[1])<1e-9 && fabs(lo[2]-hi[1]-8)<1e-9 && fabs(lo[0]-hi[2]-8)<1e-9;
        }
        fprintf(stderr,"  winding ribbon packing %d: %s, %s\n",test,stats.reason,okay?"ok":"FAIL");
        fail+=!okay;Arena_dispose(&arena);
    }
    return fail;
}

static int anp_selftest(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1),NC=5};int failures=anp_ribbon_selftest();
    for(int test=0;test<25;test++){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.continuity_ready=1;
        AsmChart charts[NC]={{0}};AsmRelation rels[3]={{0}};AsmCorr corr[3*N]={{0}};
        float xyz[NC][3*NV],uv[NC][2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
        int align=test>=16 && test<21;double common=test==4 || test==18?1e9:0,sign=test==5 || test==17?-1:1;
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
            for(int c=0;c<NC;c++){
                xyz[c][3*v]=(float)(10*c+x);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=(float)(c==4?20:0);
                if(test==19 && c==4){xyz[c][3*v]=0;xyz[c][3*v+1]=(float)(10*c+x);xyz[c][3*v+2]=(float)(y+20);}
                uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)y;
            }
            if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
        }
        if(test==13)for(int v=0;v<NV;v++)uv[4][2*v]*=1.5f;
        if(test==14 || test==20)uv[4][2*(4*N+4)]+=3;
        for(int c=0;c<NC;c++){
            AsmChart *ch=charts+c;ch->id=c;ch->parent=-1;ch->nv=NV;ch->nf=NF;
            ch->xyz=xyz[c];ch->uv=uv[c];ch->faces=faces;ch->boundary=boundary;ch->area3d=64;
            if(sign<0)ch->flags|=ASM_CHART_MIRROR;
            if(test!=15 && (test!=1 || c!=1)){ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;}
            ch->pose_x=common+100+sign*(c==4?0:10*c);ch->pose_y=common+200;
            if(align){double x=ch->pose_x-(common+100);ch->pose_x=common+100+cos(.37)*x;ch->pose_y=common+200+sin(.37)*x;ch->pose_theta=.37;}
        }
        if(test==3)charts[3].pose_x+=1000;
        for(int i=0;i<3;i++){
            AsmRelation *r=rels+i;r->a=i;r->b=i+1;r->continuity=ASM_CONT_SOURCE;
            r->corr_first=i*N;r->corr_count=N;r->n_corr=N;
            if(i==2 && (test==2 || test==3 || test==21 || test>=23))r->flags=ASM_REL_WEAK;
            for(int y=0;y<N;y++){
                AsmCorr *w=corr+i*N+y;w->va=y*N+N-1;w->vb=y*N;w->valid=3;
                w->gap_a[0]=w->gap_b[0]=2;
                w->trim_a=apr_test_trim(charts+i,w->va);w->trim_b=apr_test_trim(charts+i+1,w->vb);
            }
        }
        run.charts=charts;run.n_charts=NC;run.rels=rels;run.n_rels=3;run.corr=corr;run.n_corr=3*N;
        uint8_t member[NC]={1,1,1,1,0};
        AsmField *field=AsmField_new_local(arena,&run,test==11?member:NULL,NULL);
        if(!field){failures++;Arena_dispose(&arena);continue;}
        if(test==9)field->uv[0]=NAN;
        if(test==10 && AsmField_build_contacts(field)){failures++;Arena_dispose(&arena);continue;}
        size_t bytes=2*field->nv*sizeof(double);double *saved=ARENA_ALLOC(arena,bytes?bytes:sizeof(double));
        memcpy(saved,field->uv,bytes);
        AsmChart old_charts[NC];AsmRelation old_rels[3];AsmCorr old_corr[3*N];
        float old_xyz[NC][3*NV],old_uv[NC][2*NV];int32_t old_faces[3*NF];
        memcpy(old_charts,charts,sizeof charts);memcpy(old_rels,rels,sizeof rels);memcpy(old_corr,corr,sizeof corr);
        memcpy(old_xyz,xyz,sizeof xyz);memcpy(old_uv,uv,sizeof uv);memcpy(old_faces,faces,sizeof faces);
        double gap=test==7?0:test==8?NAN:test==12?DBL_MAX:10;
        size_t owners[NC]={0,0,0,3,4};
        if(test==23)owners[3]=NC;
        if(test==24){owners[0]=1;owners[1]=0;}
        AnpStats stats;int kept=anp_pack_layout(&run,field,gap,ves_clock_sec()+(test==6?-1:5),&stats,NULL,align,NULL,test>=21?owners:NULL);
        int expect=test<6 || test==13 || (test>=16 && test<20) || test==21;
        int okay=kept==expect && !memcmp(old_rels,rels,sizeof rels) && !memcmp(old_corr,corr,sizeof corr) &&
            !memcmp(old_xyz,xyz,sizeof xyz) && !memcmp(old_uv,uv,sizeof uv) && !memcmp(old_faces,faces,sizeof faces);
        if(expect){
            okay=okay && stats.regions==(size_t)(test==3 || test==21?3:2) && stats.charts==(size_t)(test==1?4:5) &&
                stats.faces==stats.charts*NF && stats.passing_before==stats.passing_after &&
                stats.passing_before==(size_t)(test==1?1:test==3 || test==21?2:3);
            double group_shift[2]={0,0};
            for(int c=0;c<NC;c++){
                AsmChart expected=old_charts[c];size_t first=field->vertex_offset[c];
                if(first!=SIZE_MAX){
                    expected.placed_uv=field->uv+2*first;
                    if(align){
                        size_t start=2*first;double tolerance=test==18?2e-6:1e-9;
                        for(int v=0;v<NV;v++){
                            size_t at=2*(first+(size_t)v);
                            double old=hypot(saved[at]-saved[start],saved[at+1]-saved[start+1]);
                            double next=hypot(field->uv[at]-field->uv[start],field->uv[at+1]-field->uv[start+1]);
                            okay=okay && fabs(old-next)<tolerance;
                            if(test==19 && c==4)for(int d=0;d<2;d++)okay=okay &&
                                fabs((field->uv[at+(size_t)d]-field->uv[start+(size_t)d])-(saved[at+(size_t)d]-saved[start+(size_t)d]))<tolerance;
                        }
                        if(test!=19 || c<4)okay=okay && fabs(field->uv[start+2*(N-1)]-field->uv[start])<tolerance;
                        if(c<4){
                            double old=hypot(saved[start]-saved[0],saved[start+1]-saved[1]);
                            double next=hypot(field->uv[start]-field->uv[0],field->uv[start+1]-field->uv[1]);
                            okay=okay && fabs(old-next)<tolerance;
                        }
                    }else for(int d=0;d<2;d++){
                        double shift=field->uv[2*first+d]-saved[2*first+d];
                        if(!c)group_shift[d]=shift;
                        if(c<3 || (c==3 && test!=3 && test!=21))okay=okay && shift==group_shift[d];
                        for(int v=0;v<NV;v++)okay=okay && field->uv[2*(first+(size_t)v)+d]-saved[2*(first+(size_t)v)+d]==shift;
                    }
                }
                okay=okay && !memcmp(&expected,charts+c,sizeof expected);
            }
            AsmContactStats before={0},after={0};
            AsmContacts *contacts=AsmContacts_new_threads(arena,field->nv,field->faces,field->nf,saved,1e-7,1);
            if(!contacts || AsmContacts_measure(contacts,saved,NULL,NULL,NULL,&before) ||
               AsmContacts_measure(contacts,field->uv,NULL,NULL,NULL,&after))okay=0;
            else okay=okay && before.pairs>0 && after.complete && !after.pairs;
        }else okay=okay && !memcmp(saved,field->uv,bytes) && !memcmp(old_charts,charts,sizeof charts);
        const char *reason=test==6?"deadline":test==9?"nonfinite_uv":test==11?"incomplete_field":
            test==12?"nonfinite_layout":test==14 || test==20?"metric_regression":
            test==22?"required_cross_frame":test>=23?"invalid_ownership":expect?"kept":"invalid_input";
        okay=okay && !strcmp(stats.reason,reason);
        fprintf(stderr,"  native source region packing %d: %s, %zu regions, seams %zu/%zu, %s\n",
            test,stats.reason,stats.regions,stats.passing_before,stats.passing_after,okay?"ok":"FAIL");
        failures+=!okay;Arena_dispose(&arena);
    }
    return failures;
}
