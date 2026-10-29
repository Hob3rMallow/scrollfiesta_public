/* Source-connected pieces share a frame even through missing charts.
 * Unrelated frames can be translated apart without changing a source seam. */
typedef struct AnpRegion {
    size_t root,charts,faces;
    double lo[2],hi[2],shift[2];
    double origin[2],gradient[2],cosine,sine;
    double turn_lo,turn_hi,angular_gradient[2];
    size_t supported_charts;
    /* Handedness against the scroll (winding ribbon only): area-weighted
     * sign((dP/du x dP/dv) . inward), cross product in the (z,y,x) order of
     * chart xyz -- > 0 is the handedness of the team's VC3D reading canvases
     * (sheet_layers' census, checked on PHerc0139 against the exact 3-D map
     * onto their canvases).  A rotation cannot change it: a piece flattened
     * mirror-image is presented with u -> -u (mirror = 1) and its charts'
     * ASM_CHART_MIRROR toggled, so every piece shows the recto the same way
     * round and the scroll axis runs the same way down every piece. */
    double handed;
    int mirror;
} AnpRegion;
typedef struct AnpStats {
    size_t regions,charts,faces,passing_before,passing_after,mirrored;
    int32_t first_bad_chart;
    double seconds,width,height,gap;
    int kept;
    const char *reason;
} AnpStats;

static int anp_order(const void *a,const void *b)
{
    const AnpRegion *x=a,*y=b;double hx=x->hi[1]-x->lo[1],hy=y->hi[1]-y->lo[1];
    if(hx!=hy)return hx>hy?-1:1;
    double wx=x->hi[0]-x->lo[0],wy=y->hi[0]-y->lo[0];
    if(wx!=wy)return wx>wy?-1:1;
    return x->root<y->root?-1:x->root>y->root;
}

/* The region's presentation frame of a UV point: about its origin, reflected
 * (u -> -u) for a mirror-imaged region, then rotated. */
static inline void anp_frame(const AnpRegion *r,const double *uv,double out[2])
{
    double x=uv[0]-r->origin[0],y=uv[1]-r->origin[1];
    if(r->mirror)x=-x;
    out[0]=r->cosine*x-r->sine*y;out[1]=r->sine*x+r->cosine*y;
}

/* Signed handedness weight of one chart's placed faces against the inward
 * direction at the chart's centroid (see AnpRegion.handed). */
static double anp_chart_handedness(const AsmChart *ch,const double *uv,const AsmAxis *axis)
{
    double foot[3],in[3],l=0,sum=0;
    AsmAxis_project(axis,ch->centroid,NULL,NULL,NULL,foot);
    for(int d=0;d<3;d++){in[d]=foot[d]-ch->centroid[d];l+=in[d]*in[d];}
    l=sqrt(l);if(!(l>1e-9) || !isfinite(l))return 0;
    for(int d=0;d<3;d++)in[d]/=l;
    for(size_t k=0;k<ch->nf;k++){
        const int32_t *t=ch->faces+3*k;
        const double *p=uv+2*(size_t)t[0],*q=uv+2*(size_t)t[1],*s=uv+2*(size_t)t[2];
        double det=(q[0]-p[0])*(s[1]-p[1])-(s[0]-p[0])*(q[1]-p[1]),e[3],g[3],c[3];
        if(!(fabs(det)>1e-20))continue;
        for(int d=0;d<3;d++){
            e[d]=(double)ch->xyz[3*(size_t)t[1]+(size_t)d]-ch->xyz[3*(size_t)t[0]+(size_t)d];
            g[d]=(double)ch->xyz[3*(size_t)t[2]+(size_t)d]-ch->xyz[3*(size_t)t[0]+(size_t)d];
        }
        c[0]=e[1]*g[2]-e[2]*g[1];c[1]=e[2]*g[0]-e[0]*g[2];c[2]=e[0]*g[1]-e[1]*g[0];
        sum+=(det>0?1:-1)*(c[0]*in[0]+c[1]*in[1]+c[2]*in[2]);
    }
    return isfinite(sum)?sum:0;
}

static int anp_winding_order(const void *a,const void *b)
{
    const AnpRegion *x=a,*y=b;
    if(x->turn_lo!=y->turn_lo)return x->turn_lo<y->turn_lo?-1:1;
    if(x->turn_hi!=y->turn_hi)return x->turn_hi<y->turn_hi?-1:1;
    return x->root<y->root?-1:x->root>y->root;
}

/* owners, when supplied, names each chart's certified working frame by a
 * representative chart. Independent local coordinates cannot establish a
 * join by coincidence. NULL retains the incumbent seam-derived frames. */
static int anp_pack_layout(AsmRun *run,AsmField *f,double gap,double deadline,AnpStats *out,FILE *ledger,int align,const AsmAxis *axis,const size_t *owners)
{
    if(!out)return 0;
    double started=ves_clock_sec();memset(out,0,sizeof *out);out->gap=gap;out->reason="invalid_input";out->first_bad_chart=-1;
    Arena_T arena=Arena_new();size_t nc=run?run->n_charts:0;
    if(!run || !f || !nc || !f->nv || !(gap>0) || !isfinite(gap) || f->contacts || f->run!=run)goto done;
    for(size_t c=0;c<nc;c++)if(AsmChart_registered(run->charts+c)!=(f->vertex_offset[c]!=SIZE_MAX)){
        out->reason="incomplete_field";goto done;
    }
    if(deadline>0 && ves_clock_sec()>=deadline){out->reason="deadline";goto done;}
    AsmRibbonKey *keys=NULL;
    if(axis){
        keys=ARENA_CALLOC(arena,nc,sizeof *keys);
        if(AsmReadingOrder_ribbon_keys(run,axis,keys)){out->reason="winding_evidence_failed";goto done;}
        align=1;
    }
    size_t *group=ARENA_ALLOC(arena,nc*sizeof *group),*index=ARENA_ALLOC(arena,nc*sizeof *index);
    for(size_t c=0;c<nc;c++){
        group[c]=owners?owners[c]:c;index[c]=SIZE_MAX;
        if(owners && (owners[c]>=nc || owners[owners[c]]!=owners[c] ||
           (f->vertex_offset[c]!=SIZE_MAX && f->vertex_offset[owners[c]]==SIZE_MAX))){
            out->reason="invalid_ownership";goto done;
        }
    }
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *s=f->seams+i;AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);
        if(owners){
            int same=group[s->a]==group[s->b];
            if(!same && s->required && f->vertex_offset[s->a]!=SIZE_MAX && f->vertex_offset[s->b]!=SIZE_MAX){
                out->reason="required_cross_frame";goto done;
            }
            out->passing_before+=same && m.pass;
        }else{
            out->passing_before+=m.pass;
            if(s->required || m.pass){size_t a=apr_root(group,s->a),b=apr_root(group,s->b);group[a>b?a:b]=a<b?a:b;}
        }
    }
    AnpRegion *regions=ARENA_CALLOC(arena,nc,sizeof *regions);
    for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        size_t root=apr_root(group,c),at=index[root];
        if(at==SIZE_MAX){
            at=index[root]=out->regions++;regions[at].root=root;regions[at].cosine=1;
            regions[at].turn_lo=DBL_MAX;regions[at].turn_hi=-DBL_MAX;
            for(int d=0;d<2;d++){regions[at].lo[d]=DBL_MAX;regions[at].hi[d]=-DBL_MAX;}
        }
        AnpRegion *r=regions+at;r->charts++;r->faces+=run->charts[c].nf;
        if(keys){
            r->turn_lo=fmin(r->turn_lo,keys[c].turn_lo);r->turn_hi=fmax(r->turn_hi,keys[c].turn_hi);
            r->supported_charts+=(size_t)keys[c].source_supported;
            for(int d=0;d<2;d++)r->angular_gradient[d]+=keys[c].angular_gradient[d];
            r->handed+=anp_chart_handedness(run->charts+c,f->uv+2*f->vertex_offset[c],axis);
        }
        out->charts++;out->faces+=run->charts[c].nf;
        for(size_t v=0;v<run->charts[c].nv;v++)for(int d=0;d<2;d++){
            double x=f->uv[2*(f->vertex_offset[c]+v)+d];
            if(!isfinite(x)){out->reason="nonfinite_uv";goto done;}
            r->lo[d]=fmin(r->lo[d],x);r->hi[d]=fmax(r->hi[d],x);
        }
        /* Integrate the source scan-Z gradient over represented UV faces.
         * This chooses a presentation frame, never an elastic deformation. */
        if(align){
            const AsmChart *ch=run->charts+c;size_t first=f->vertex_offset[c];
            for(size_t k=0;k<ch->nf;k++){
                const int32_t *t=ch->faces+3*k;
                const double *p=f->uv+2*(first+(size_t)t[0]),*q=f->uv+2*(first+(size_t)t[1]),*s=f->uv+2*(first+(size_t)t[2]);
                double u=q[0]-p[0],v=q[1]-p[1],a=s[0]-p[0],b=s[1]-p[1],det=u*b-v*a;
                double dz=(double)ch->xyz[3*(size_t)t[1]]-ch->xyz[3*(size_t)t[0]],ez=(double)ch->xyz[3*(size_t)t[2]]-ch->xyz[3*(size_t)t[0]];
                if(fabs(det)>1e-20){double sign=det>0?1:-1;r->gradient[0]+=.5*sign*(dz*b-ez*v);r->gradient[1]+=.5*sign*(ez*u-dz*a);}
            }
        }
    }
    if(align){
        for(size_t i=0;i<out->regions;i++){
            AnpRegion *r=regions+i;
            if(keys && r->handed<0){        /* present the mirror image's reflection: u -> -u */
                r->mirror=1;out->mirrored++;
                r->gradient[0]=-r->gradient[0];r->angular_gradient[0]=-r->angular_gradient[0];
            }
            double length=hypot(r->gradient[0],r->gradient[1]);
            if(!isfinite(length)){out->reason="nonfinite_axis_gradient";goto done;}
            if(length>1e-20){
                r->cosine=r->gradient[1]/length;r->sine=r->gradient[0]/length;
                /* Choose the equivalent frame with at most a quarter turn.
                 * Scan-Z alignment does not determine reading direction. */
                if(r->cosine<0){r->cosine=-r->cosine;r->sine=-r->sine;}
            }
            if(keys){
                double phase_length=hypot(r->angular_gradient[0],r->angular_gradient[1]);
                if(phase_length>1e-20){
                    r->cosine=r->angular_gradient[0]/phase_length;
                    r->sine=-r->angular_gradient[1]/phase_length;
                }
            }
            for(int d=0;d<2;d++){r->origin[d]=r->lo[d];r->lo[d]=DBL_MAX;r->hi[d]=-DBL_MAX;}
        }
        for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX){
            AnpRegion *r=regions+index[apr_root(group,c)];size_t first=f->vertex_offset[c];
            for(size_t v=0;v<run->charts[c].nv;v++){
                double point[2];anp_frame(r,f->uv+2*(first+v),point);
                for(int d=0;d<2;d++){r->lo[d]=fmin(r->lo[d],point[d]);r->hi[d]=fmax(r->hi[d],point[d]);}
            }
        }
    }
    qsort(regions,out->regions,sizeof *regions,keys?anp_winding_order:anp_order);
    double area=0,wide=0;
    for(size_t i=0;i<out->regions;i++){
        AnpRegion *r=regions+i;index[r->root]=i;
        double w=r->hi[0]-r->lo[0],h=r->hi[1]-r->lo[1];area+=(w+gap)*(h+gap);wide=fmax(wide,w);
        if(!isfinite(area)){out->reason="nonfinite_layout";goto done;}
    }
    double limit=fmax(wide,sqrt(area)),x=0,y=0,row_height=0;
    for(size_t i=0;i<out->regions;i++){
        AnpRegion *r=regions+i;double w=r->hi[0]-r->lo[0],h=r->hi[1]-r->lo[1];
        if(!keys && x>0 && x+w>limit){y+=row_height+gap;x=0;row_height=0;}
        r->shift[0]=x-r->lo[0];r->shift[1]=y-r->lo[1];
        out->width=fmax(out->width,x+w);out->height=fmax(out->height,y+h);
        x+=w+gap;row_height=fmax(row_height,h);
        if(!isfinite(x) || !isfinite(y) || !isfinite(out->width) || !isfinite(out->height)){out->reason="nonfinite_layout";goto done;}
    }
    double *candidate=ARENA_ALLOC(arena,2*f->nv*sizeof(double));
    for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        const AnpRegion *r=regions+index[apr_root(group,c)];size_t first=f->vertex_offset[c];
        for(size_t v=0;v<run->charts[c].nv;v++){
            double point[2];anp_frame(r,f->uv+2*(first+v),point);
            candidate[2*(first+v)]=point[0]+r->shift[0];
            candidate[2*(first+v)+1]=point[1]+r->shift[1];
        }
    }
    /* Certify every original face before publishing rigidly moved coordinates,
     * including inherited imperfect charts and parity-encoded source UVs. */
    for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        if(deadline>0 && ves_clock_sec()>=deadline){out->reason="deadline";goto done;}
        Arena_Mark mark=Arena_save(arena);AsmChart chart=run->charts[c];
        chart.placed_uv=f->uv+2*f->vertex_offset[c];
        AsmMetricFace *metric=ARENA_ALLOC(arena,chart.nf*sizeof *metric);AsmMetricStats old={0},next={0};
        if(AsmMetric_prepare(&chart,metric)){out->reason="source_metric_failed";out->first_bad_chart=(int32_t)c;goto done;}
        /* a reflected region is certified under the parity it will carry */
        AsmChart moved=chart;if(regions[index[apr_root(group,c)]].mirror)moved.flags^=ASM_CHART_MIRROR;
        AsmMetric_measure_encoded_uv64(&chart,metric,chart.placed_uv,&old);
        AsmMetric_measure_encoded_uv64(&moved,metric,candidate+2*f->vertex_offset[c],&next);
        if(!AsmMetric_preserved(&old,&next) ||
           !AsmMetric_encoded_faces_preserved_parity(&chart,metric,candidate+2*f->vertex_offset[c],
                                                     (moved.flags&ASM_CHART_MIRROR)!=0)){
            out->reason="metric_regression";out->first_bad_chart=(int32_t)c;goto done;
        }
        Arena_restore(arena,mark);
    }
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *seam=f->seams+i;
        int same=!owners || group[seam->a]==group[seam->b];
        AsmSeamMeasure old,next;AsmField_seam(f,i,f->uv,&old);AsmField_seam(f,i,candidate,&next);
        out->passing_after+=same && next.pass;
        if((same && old.pass && !next.pass) || (seam->required &&
           (old.invalid!=next.invalid || next.energy>old.energy+1e-8*fmax(1,old.energy)))){
            out->reason="seam_regression";goto done;
        }
    }
    if(deadline>0 && ves_clock_sec()>=deadline){out->reason="deadline";goto done;}
    if(ledger){
        if(fputs("chart,protected_region,du,dv,origin_u,origin_v,cosine,sine,sequence,approx_turn_lo,approx_turn_hi,source_supported_charts,mirror\n",ledger)<0){out->reason="ledger_failed";goto done;}
        for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX){
            const AnpRegion *r=regions+index[apr_root(group,c)];
            if(fprintf(ledger,"%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu,%.17g,%.17g,%zu,%d\n",c,r->root,r->shift[0],r->shift[1],r->origin[0],r->origin[1],r->cosine,r->sine,index[r->root],keys?r->turn_lo:0,keys?r->turn_hi:0,r->supported_charts,r->mirror)<0){out->reason="ledger_failed";goto done;}
        }
        if(fflush(ledger)){out->reason="ledger_failed";goto done;}
    }
    /* The complete field and its owner were checked above. Publishing the
     * rigid motion and its chart pointers has no remaining failure path. */
    memcpy(f->uv,candidate,2*f->nv*sizeof(double));
    for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        run->charts[c].placed_uv=f->uv+2*f->vertex_offset[c];
        if(regions[index[apr_root(group,c)]].mirror)run->charts[c].flags^=ASM_CHART_MIRROR;
    }
    out->kept=1;out->reason="kept";
done:
    out->seconds=ves_clock_sec()-started;Arena_dispose(&arena);return out->kept;
}
