/* Read-only source-to-atlas index. Missing-chart pins are averages of their
 * represented boundary observations, not invented placements or hole fills. */
static int asv_write(const AsmRun *run,const char *dir)
{
    Arena_T arena=Arena_new();size_t nc=run->n_charts;int rc=-1;
    size_t *group=ARENA_ALLOC(arena,nc*sizeof *group);
    double *pin=ARENA_CALLOC(arena,3*nc,sizeof *pin);
    for(size_t c=0;c<nc;c++)group[c]=c;
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;if(!apa_required_source(r))continue;
        size_t a=apr_root(group,r->a),b=apr_root(group,r->b);group[a>b?a:b]=a<b?a:b;
        int pa=AsmChart_registered(run->charts+r->a),pb=AsmChart_registered(run->charts+r->b);
        if(pa==pb)continue;
        size_t missing=(size_t)(pa?r->b:r->a);const AsmChart *ch=run->charts+(pa?r->a:r->b);
        for(int k=0;k<r->corr_count;k++){
            const AsmCorr *w=run->corr+r->corr_first+k;int32_t v=pa?w->va:w->vb;
            if(v<0 || (size_t)v>=ch->nv)continue;double p[2];AsmChart_point(ch,(size_t)v,p);
            pin[3*missing]+=p[0];pin[3*missing+1]+=p[1];pin[3*missing+2]++;
        }
    }
    char path[2048];if(snprintf(path,sizeof path,"%s/sheet_index.json",dir)>=(int)sizeof path)goto done;
    ves_mkdir(dir);FILE *fp=fopen(path,"wb");if(!fp)goto done;
    fprintf(fp,"{\"schema\":\"native-sheet-review-v1\",\"source_modified\":false,\"pin_rule\":\"mean represented required-boundary observations; not a proposed placement\",\"charts\":[\n");
    for(size_t c=0;c<nc;c++){
        const AsmChart *ch=run->charts+c;int placed=AsmChart_registered(ch);
        double low[2]={DBL_MAX,DBL_MAX},high[2]={-DBL_MAX,-DBL_MAX},sum[5]={0},gradient[2]={0};
        double point[16][2]={{0}},extreme[16];for(int d=0;d<16;d++)extreme[d]=-DBL_MAX;
        if(placed)for(size_t v=0;v<ch->nv;v++){
            double p[2];AsmChart_point(ch,v,p);
            for(int d=0;d<2;d++){low[d]=fmin(low[d],p[d]);high[d]=fmax(high[d],p[d]);sum[d]+=p[d];}
            sum[2]+=p[0]*p[0];sum[3]+=p[0]*p[1];sum[4]+=p[1]*p[1];
            for(int d=0;d<16;d++){
                double angle=6.283185307179586*(double)d/16,value=p[0]*cos(angle)+p[1]*sin(angle);
                if(value>extreme[d]){extreme[d]=value;point[d][0]=p[0];point[d][1]=p[1];}
            }
        }
        if(placed)for(size_t f=0;f<ch->nf;f++){
            const int32_t *t=ch->faces+3*f;double p[3][2];for(int k=0;k<3;k++)AsmChart_point(ch,(size_t)t[k],p[k]);
            double u=p[1][0]-p[0][0],v=p[1][1]-p[0][1],a=p[2][0]-p[0][0],b=p[2][1]-p[0][1];
            double det=u*b-v*a,dz=(double)ch->xyz[3*(size_t)t[1]]-ch->xyz[3*(size_t)t[0]],ez=(double)ch->xyz[3*(size_t)t[2]]-ch->xyz[3*(size_t)t[0]];
            if(fabs(det)>1e-20){double sign=det>0?1:-1;gradient[0]+=.5*sign*(dz*b-ez*v);gradient[1]+=.5*sign*(ez*u-dz*a);}
        }
        const char *cube=ch->cube>=0 && (size_t)ch->cube<run->n_cubes?run->pile[ch->cube].cube_id:"";
        fprintf(fp,"%s{\"id\":%zu,\"region\":%zu,\"cube\":\"%s\",\"flags\":%u,\"placed\":%s,\"nv\":%zu,\"faces\":%zu,\"area\":%.17g,\"xyz_lo\":[%.9g,%.9g,%.9g],\"xyz_hi\":[%.9g,%.9g,%.9g],\"uv_box\":",
            c?",\n":"",c,apr_root(group,c),cube,ch->flags,placed?"true":"false",ch->nv,ch->nf,ch->area3d,
            (double)ch->bbox_lo[0],(double)ch->bbox_lo[1],(double)ch->bbox_lo[2],(double)ch->bbox_hi[0],(double)ch->bbox_hi[1],(double)ch->bbox_hi[2]);
        if(placed)fprintf(fp,"[%.17g,%.17g,%.17g,%.17g]",low[0],low[1],high[0],high[1]);else fputs("null",fp);
        fprintf(fp,",\"moments\":[%.17g,%.17g,%.17g,%.17g,%.17g],\"z_gradient_integral\":[%.17g,%.17g],\"hull\":[",sum[0],sum[1],sum[2],sum[3],sum[4],gradient[0],gradient[1]);
        if(placed)for(int d=0;d<16;d++)fprintf(fp,"%s[%.17g,%.17g]",d?",":"",point[d][0],point[d][1]);
        fputs("],\"missing_boundary_pin\":",fp);
        if(pin[3*c+2]>0)fprintf(fp,"[%.17g,%.17g]",pin[3*c]/pin[3*c+2],pin[3*c+1]/pin[3*c+2]);else fputs("null",fp);
        fputs("}",fp);
    }
    fputs("\n]}\n",fp);int okay=!ferror(fp);if(fclose(fp))okay=0;if(okay)rc=0;
done:
    Arena_dispose(&arena);return rc;
}

static int asv_case(const char *input,const char *dir,const char *axis_path)
{
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;int rc=-1;
    AsmAxis *axis=axis_path?AsmAxis_load(arena,axis_path):NULL;
    if(axis_path && !axis)goto done;
    fprintf(stderr,"[sheet review] loading checkpoint\n");
    if(!AsmStore_read_run(input,&run)){
        fprintf(stderr,"[sheet review] loaded %zu charts; writing source index\n",run.n_charts);
        rc=asv_write(&run,dir);
        if(!rc && axis)rc=AsmReadingOrder_write(&run,axis,dir);
    }
done:
    Arena_dispose(&arena);return rc;
}

int AsmRepair_sheet_review_case(const char *input,const char *dir)
{return asv_case(input,dir,NULL);}

int AsmRepair_reading_review_case(const char *input,const char *dir,const char *axis_path)
{return asv_case(input,dir,axis_path);}

/* Pieces for a presentation pack: charts joined by seams between two placed
 * charts, passing or required. A seam to an unplaced chart ties nothing:
 * separately placed pieces are independent frames, and grouping through a
 * missing chart glued them at an arbitrary offset (PHerc0139 10^3 r51 review:
 * a 16-chart piece 8,500 vox below the sheet, tied to the largest piece by two
 * required seams of one unplaced chart). Each chart is owned by the smallest
 * chart index of its piece, as anp_pack_layout's owners require. */
static size_t *asv_piece_owners(Arena_T arena,const AsmRun *run,const AsmField *field)
{
    size_t *owner=ARENA_ALLOC(arena,(run->n_charts?run->n_charts:1)*sizeof *owner);
    for(size_t c=0;c<run->n_charts;c++)owner[c]=c;
    for(size_t i=0;i<field->n_seams;i++){
        const AsmFieldSeam *s=field->seams+i;
        if(field->vertex_offset[s->a]==SIZE_MAX || field->vertex_offset[s->b]==SIZE_MAX)continue;
        if(!s->required){AsmSeamMeasure m;AsmField_seam(field,i,field->uv,&m);if(!m.pass)continue;}
        size_t a=apr_root(owner,(size_t)s->a),b=apr_root(owner,(size_t)s->b);owner[a>b?a:b]=a<b?a:b;
    }
    for(size_t c=0;c<run->n_charts;c++)owner[c]=apr_root(owner,c);
    return owner;
}

/* After anp_pack_layout every placed chart's placed_uv points into the pack
 * field's coordinates. Each gets its own copy in `arena`, so the field -- a
 * full copy of the run's geometry, ~38 GB at PHerc0139 21^3 -- can be released
 * before the audit builds its own. Both stayed live on 2026-09-29, and the
 * 21^3 ribbon's audit helped exhaust the 272 GB commit limit. */
static void asv_own_placed_uv(Arena_T arena,AsmRun *run,const AsmField *f)
{
    for(size_t c=0;c<run->n_charts;c++){
        AsmChart *ch=run->charts+c;
        if(f->vertex_offset[c]==SIZE_MAX || !ch->placed_uv)continue;
        double *uv=ARENA_ALLOC(arena,2*ch->nv*sizeof *uv);
        memcpy(uv,ch->placed_uv,2*ch->nv*sizeof *uv);ch->placed_uv=uv;
    }
}

int AsmRepair_winding_ribbon_case(const char *input,const char *axis_path,const char *dir)
{
    double started=ves_clock_sec();Arena_T arena=Arena_new(),field_arena=NULL;AsmRun run={0};run.arena=arena;int rc=-1;
    AsmAxis *axis=AsmAxis_load(arena,axis_path);FILE *file=NULL;char path[2048];
    if(!axis || AsmStore_read_run(input,&run))goto done;
    AsmField_threads(24);field_arena=Arena_new();AsmField *field=AsmField_new_local(field_arena,&run,NULL,NULL);
    if(!field)goto done;
    ves_mkdir(dir);snprintf(path,sizeof path,"%s/region_pack.csv",dir);file=fopen(path,"wb");if(!file)goto done;
    AnpStats packed;int okay=anp_pack_layout(&run,field,ASM_READING_RIBBON_GAP,0,&packed,file,1,axis,asv_piece_owners(arena,&run,field));
    if(fclose(file))okay=0;file=NULL;
    snprintf(path,sizeof path,"%s/ribbon_layout.json",dir);file=fopen(path,"wb");if(!file)goto done;
    if(fprintf(file,"{\"schema\":\"native-approximate-winding-ribbon-v1\",\"kept\":%s,\"reason\":\"%s\","
        "\"ordering\":\"approximate_winding\",\"identity_certified\":false,\"connected_strip\":false,"
        "\"scope\":\"One horizontal sequence of intact protected pieces, ordered by source phase with an approximate radius gauge. Gaps separate disconnected pieces; no material is added.\","
        "\"regions\":%zu,\"charts\":%zu,\"faces\":%zu,\"pitch\":%.17g,\"gap\":%.17g,\"width\":%.17g,\"height\":%.17g,"
        "\"passing_before\":%zu,\"passing_after\":%zu,\"mirrored_regions\":%zu,"
        "\"handedness\":\"every region presented with the team's reading-canvas handedness: (dP/du x dP/dv) . inward > 0, (z,y,x) order; mirror-imaged regions reflected u -> -u with ASM_CHART_MIRROR toggled\","
        "\"pack_seconds\":%.17g}\n",
        packed.kept?"true":"false",packed.reason,packed.regions,packed.charts,packed.faces,run.layer_d_global,
        packed.gap,packed.width,packed.height,packed.passing_before,packed.passing_after,packed.mirrored,
        packed.seconds)<0)okay=0;
    if(fclose(file))okay=0;file=NULL;
    if(!okay){
        fprintf(stderr,"[winding ribbon] refused: %s, first bad chart %d\n",packed.reason,(int)packed.first_bad_chart);
        if(packed.first_bad_chart>=0 && (size_t)packed.first_bad_chart<run.n_charts){
            const AsmChart *b=run.charts+packed.first_bad_chart;
            fprintf(stderr,"[winding ribbon] chart %d: cube %d comp %d parent %d flags 0x%x nv %zu nf %zu area %.1f sigma [%.6f, %.6f] stress %.4f\n",
                    b->id,b->cube,b->comp,b->parent,b->flags,b->nv,b->nf,b->area3d,b->sigma_lo,b->sigma_hi,b->stress_frac);
        }
        goto done;
    }
    fprintf(stderr,"[winding ribbon] %zu intact pieces (%zu reflected to the reading handedness), %.3f x %.3f vox; %.3f seconds\n",
            packed.regions,packed.mirrored,packed.width,packed.height,ves_clock_sec()-started);
    asv_own_placed_uv(arena,&run,field);Arena_dispose(&field_arena);field=NULL;
    AsmAuditStats audit={0};
    if(AsmAudit_checkpoint(&run,dir,&audit) || !audit.complete || !audit.source_preserved || !audit.contacts_complete ||
       audit.invalid_faces || audit.strict_bad_faces || audit.bad_metric_charts || audit.overlapping_pairs || audit.unresolved_present_repair_seams)goto done;
    snprintf(path,sizeof path,"%s/proposal.asr",dir);
    if(AsmStore_write_run(path,&run) || asv_write(&run,dir))goto done;
    rc=0;
done:
    if(file)fclose(file);if(field_arena)Arena_dispose(&field_arena);Arena_dispose(&arena);return rc;
}

static int asv_selftest(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
    AsmChart chart={0};float xyz[]={0,0,0,1,0,0,0,1,0},uv[]={0,0,1,0,0,1};int32_t faces[]={0,1,2};
    chart.id=0;chart.cube=-1;chart.nv=3;chart.nf=1;chart.xyz=xyz;chart.uv=uv;chart.faces=faces;
    chart.placed=1;chart.placement_state=ASM_PLACE_ROOT;chart.area3d=.5;chart.bbox_hi[0]=chart.bbox_hi[1]=1;
    run.charts=&chart;run.n_charts=1;AsmChart saved=chart;
    ves_mkdir("output");int fail=asv_write(&run,"output/sheet_review_selftest") || memcmp(&chart,&saved,sizeof chart);
    FILE *fp=fopen("output/sheet_review_selftest/sheet_index.json","rb");char buffer[4096]={0};
    if(!fp)fail=1;else{size_t n=fread(buffer,1,sizeof buffer-1,fp);buffer[n]=0;fclose(fp);
        fail|=!strstr(buffer,"\"z_gradient_integral\":[0.5,0]") || !strstr(buffer,"\"uv_box\":[0,0,1,1]");}
    {   /* Presentation pieces join placed charts only: with the grid's middle
         * column unplaced, the left and right columns are two pieces, although
         * required seams run through the middle column. */
        Arena_T grid=Arena_new();AsmRun g={0};int bad=1;
        if(!ags_grid_fixture(grid,&g)){
            for(size_t c=0;c<g.n_charts;c++){
                AsmChart *ch=g.charts+c;int use=c%3!=1;
                ch->placed=use;ch->placement_state=use?ASM_PLACE_ROOT:ASM_PLACE_NONE;ch->placed_uv=NULL;
                if(!use)continue;
                ch->placed_uv=ARENA_ALLOC(grid,2*ch->nv*sizeof(double));
                for(size_t v=0;v<ch->nv;v++){ch->placed_uv[2*v]=ch->xyz[3*v];ch->placed_uv[2*v+1]=ch->xyz[3*v+1];}
            }
            AsmField *gf=AsmField_new_local(grid,&g,NULL,NULL);
            if(gf){
                size_t *owner=asv_piece_owners(grid,&g,gf);AnpStats st;
                bad=owner[0]!=owner[3] || owner[3]!=owner[6] || owner[2]!=owner[5] || owner[5]!=owner[8] || owner[0]==owner[2] || owner[1]!=1;
                bad|=!anp_pack_layout(&g,gf,8,0,&st,NULL,0,NULL,owner) || st.regions!=2;
            }
        }
        fail|=bad;fprintf(stderr,"  native presentation pieces join placed charts only: %s\n",bad?"FAIL":"ok");Arena_dispose(&grid);
    }
    {   /* The ribbon releases its pack field before the audit: after the pack
         * the placed charts' coordinates alias the field; asv_own_placed_uv
         * copies them out, and they survive the field's arena being disposed. */
        Arena_T grid=Arena_new(),fa=Arena_new();AsmRun g={0};int bad=1;
        if(!ags_grid_fixture(grid,&g)){
            size_t total=0;
            for(size_t c=0;c<g.n_charts;c++){
                AsmChart *ch=g.charts+c;ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;total+=2*ch->nv;
                ch->placed_uv=ARENA_ALLOC(grid,2*ch->nv*sizeof(double));
                for(size_t v=0;v<ch->nv;v++){ch->placed_uv[2*v]=ch->xyz[3*v];ch->placed_uv[2*v+1]=ch->xyz[3*v+1];}
            }
            AsmField *gf=AsmField_new_local(fa,&g,NULL,NULL);AnpStats st;
            if(gf && anp_pack_layout(&g,gf,8,0,&st,NULL,0,NULL,asv_piece_owners(grid,&g,gf))){
                uintptr_t lo=(uintptr_t)gf->uv,hi=(uintptr_t)(gf->uv+2*gf->nv);
                double *snap=ARENA_ALLOC(grid,(total?total:1)*sizeof *snap);size_t k=0;int before=0,after=0,same=1;
                for(size_t c=0;c<g.n_charts;c++){
                    const AsmChart *ch=g.charts+c;uintptr_t p=(uintptr_t)ch->placed_uv;before+=p>=lo && p<hi;
                    memcpy(snap+k,ch->placed_uv,2*ch->nv*sizeof *snap);k+=2*ch->nv;
                }
                asv_own_placed_uv(grid,&g,gf);
                for(size_t c=0;c<g.n_charts;c++){uintptr_t p=(uintptr_t)g.charts[c].placed_uv;after+=p>=lo && p<hi;}
                Arena_dispose(&fa);k=0;
                for(size_t c=0;c<g.n_charts;c++){
                    const AsmChart *ch=g.charts+c;same&=!memcmp(snap+k,ch->placed_uv,2*ch->nv*sizeof *snap);k+=2*ch->nv;
                }
                bad=before!=(int)g.n_charts || after || !same;
            }
        }
        if(fa)Arena_dispose(&fa);
        fail|=bad;fprintf(stderr,"  winding ribbon releases its pack field before the audit: %s\n",bad?"FAIL":"ok");Arena_dispose(&grid);
    }
    {   /* Winding ribbon handedness: two flat 10 x 10 charts facing a straight
         * axis (z through y = x = 0) from x = 50; A is placed as flattened, B is
         * placed mirror-imaged (ASM_CHART_MIRROR, u -> -u).  The ribbon must
         * present both with the team's reading handedness -- B reflected and its
         * mirror bit cleared, A untouched -- with the scroll axis running the
         * same way down both. */
        Arena_T hand=Arena_new();AsmRun g={0};g.arena=hand;g.layer_d_global=10;int bad=0;
        static float hx[2][12]={{0,0,50, 0,10,50, 10,10,50, 10,0,50},{20,0,50, 20,10,50, 30,10,50, 30,0,50}};
        static float huv[8]={0,0, 10,0, 10,10, 0,10};
        static int32_t hf[6]={0,1,2, 0,2,3};
        double point[3]={0,0,0},dir[3]={1,0,0};
        AsmAxis *ax=AsmAxis_line(hand,point,dir,-1000,1000);
        AsmChart *ch=ARENA_CALLOC(hand,2,sizeof *ch);
        for(int k=0;k<2;k++){
            AsmChart *c=ch+k;c->id=k;c->cube=-1;c->parent=-1;c->component=-1;c->nv=4;c->nf=2;
            c->xyz=hx[k];c->uv=huv;c->faces=hf;c->placed=1;c->placement_state=ASM_PLACE_ROOT;c->area3d=100;
            c->placed_uv=ARENA_ALLOC(hand,8*sizeof(double));
            for(int v=0;v<4;v++){
                double u=huv[2*v],w=huv[2*v+1];
                c->placed_uv[2*v]=k?100-u:u;c->placed_uv[2*v+1]=w+(k?20:0);
            }
            if(k)c->flags=ASM_CHART_MIRROR;
            for(int d=0;d<3;d++){
                c->centroid[d]=.25*((double)hx[k][d]+hx[k][3+d]+hx[k][6+d]+hx[k][9+d]);
                c->bbox_lo[d]=fminf(fminf(hx[k][d],hx[k][3+d]),fminf(hx[k][6+d],hx[k][9+d]));
                c->bbox_hi[d]=fmaxf(fmaxf(hx[k][d],hx[k][3+d]),fmaxf(hx[k][6+d],hx[k][9+d]));
            }
        }
        g.charts=ch;g.n_charts=2;
        double before[2]={anp_chart_handedness(ch,ch[0].placed_uv,ax),anp_chart_handedness(ch+1,ch[1].placed_uv,ax)};
        AsmField *hf2=ax?AsmField_new_local(hand,&g,NULL,NULL):NULL;AnpStats st;memset(&st,0,sizeof st);
        if(!hf2 || !(before[0]>0) || !(before[1]<0))bad=1;
        else if(!anp_pack_layout(&g,hf2,8,0,&st,NULL,1,ax,asv_piece_owners(hand,&g,hf2)))bad=1;
        else{
            double after[2]={anp_chart_handedness(ch,ch[0].placed_uv,ax),anp_chart_handedness(ch+1,ch[1].placed_uv,ax)};
            /* z runs the same way down both: the v step from vertex 0 to vertex 3 (same y, z + 10) */
            double dz0=ch[0].placed_uv[2*3+1]-ch[0].placed_uv[1],dz1=ch[1].placed_uv[2*3+1]-ch[1].placed_uv[1];
            bad=!st.kept || st.mirrored!=1 || !(after[0]>0) || !(after[1]>0) ||
                (ch[0].flags&ASM_CHART_MIRROR) || (ch[1].flags&ASM_CHART_MIRROR) || (dz0>0)!=(dz1>0);
        }
        fail|=bad;fprintf(stderr,"  winding ribbon presents every piece with the reading handedness: %s (reason %s)\n",
                          bad?"FAIL":"ok",st.reason?st.reason:"-");
        Arena_dispose(&hand);
    }
    fprintf(stderr,"  native source-linked sheet review: %s\n",fail?"FAIL":"ok");Arena_dispose(&arena);return fail;
}
