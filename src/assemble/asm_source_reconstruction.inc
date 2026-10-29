/* Source-only boundary reconstruction, corresponding to boundary_hypothesis
 * and internal_boundary_hypothesis.
 * Cube-wide nearest neighbours can mask a continuation with a third chart.
 * Search each original chart pair on its common cube plane, or its facing
 * boundary fronts inside one cube. UV fit is never a selection criterion;
 * tangent frames are frozen before placement repair. */
static int asr_domain(const AsmRun *run,const AsmChart *a,const AsmChart *b,int *axis,double *plane)
{
    if(!run->pile || !(run->cube_size>0) || a->cube<0 || b->cube<0 ||
       (size_t)a->cube>=run->n_cubes || (size_t)b->cube>=run->n_cubes)return 0;
    const MeshPileEntry *pa=run->pile+a->cube,*pb=run->pile+b->cube;
    if(!pa->has_id || !pb->has_id)return 0;
    double x[3]={(double)pa->oz,(double)pa->oy,(double)pa->ox};
    double y[3]={(double)pb->oz,(double)pb->oy,(double)pb->ox};int found=-1;
    for(int d=0;d<3;d++)if(x[d]!=y[d]){
        if(found>=0 || fabs(fabs(x[d]-y[d])-run->cube_size)>1e-8)return 0;
        found=d;
    }
    *axis=found;*plane=found<0?0:fmax(x[found],y[found]);return 1;
}
static void asr_normal(const ActFrame *f,double n[3])
{for(int q=0;q<3;q++)n[q]=f->map[0][(q+1)%3]*f->map[1][(q+2)%3]-f->map[0][(q+2)%3]*f->map[1][(q+1)%3];}

typedef struct AsrFront {double normal[3],outward[3];size_t degree;int valid;} AsrFront;
typedef struct AsrEdge {int32_t lo,hi,a,b;double outward[3];} AsrEdge;
static int asr_edge_order(const void *a,const void *b)
{
    const AsrEdge *x=a,*y=b;
    if(x->lo!=y->lo)return x->lo<y->lo?-1:1;
    return x->hi<y->hi?-1:x->hi>y->hi;
}

/* Original oriented boundary edges determine the outward tangent direction.
 * Equally directed fronts of nearby parallel layers or overlapping copies
 * are not continuations. UV reflection must not reverse this physical test. */
static AsrFront *asr_front_geometry(Arena_T arena,const AsmChart *c)
{
    if(!c->xyz || !c->faces || !c->nv || !c->nf || c->nv>INT32_MAX ||
       c->nv>SIZE_MAX/sizeof(AsrFront) || c->nf>SIZE_MAX/(3*sizeof(AsrEdge)))return NULL;
    AsrFront *front=ARENA_CALLOC(arena,c->nv,sizeof *front);
    AsrEdge *edges=ARENA_ALLOC(arena,3*c->nf*sizeof *edges);
    for(size_t t=0;t<c->nf;t++){
        const int32_t *f=c->faces+3*t;double e[2][3],normal[3];
        for(int k=0;k<3;k++)if(f[k]<0 || (size_t)f[k]>=c->nv)return NULL;
        for(int k=0;k<2;k++)for(int d=0;d<3;d++)e[k][d]=(double)c->xyz[3*(size_t)f[k+1]+d]-c->xyz[3*(size_t)f[0]+d];
        for(int d=0;d<3;d++)normal[d]=e[0][(d+1)%3]*e[1][(d+2)%3]-e[0][(d+2)%3]*e[1][(d+1)%3];
        double length=hypot(hypot(normal[0],normal[1]),normal[2]);
        if(!(length>0) || !isfinite(length))return NULL;
        for(int k=0;k<3;k++)for(int d=0;d<3;d++)front[f[k]].normal[d]+=normal[d];
        for(int k=0;k<3;k++){
            int32_t a=f[k],b=f[(k+1)%3];AsrEdge *edge=edges+3*t+k;double direction[3];
            edge->a=a;edge->b=b;edge->lo=a<b?a:b;edge->hi=a<b?b:a;
            for(int d=0;d<3;d++)direction[d]=(double)c->xyz[3*(size_t)b+d]-c->xyz[3*(size_t)a+d];
            for(int d=0;d<3;d++)edge->outward[d]=(direction[(d+1)%3]*normal[(d+2)%3]-direction[(d+2)%3]*normal[(d+1)%3])/length;
        }
    }
    qsort(edges,3*c->nf,sizeof *edges,asr_edge_order);
    for(size_t first=0;first<3*c->nf;){
        size_t last=first+1;while(last<3*c->nf && edges[last].lo==edges[first].lo && edges[last].hi==edges[first].hi)last++;
        if(last-first>2)return NULL;
        if(last-first==1)for(int k=0;k<2;k++){
            AsrFront *point=front+(k?edges[first].b:edges[first].a);point->degree++;
            for(int d=0;d<3;d++)point->outward[d]+=edges[first].outward[d];
        }
        first=last;
    }
    for(size_t v=0;v<c->nv;v++){
        AsrFront *point=front+v;
        double nl=hypot(hypot(point->normal[0],point->normal[1]),point->normal[2]);
        double ol=hypot(hypot(point->outward[0],point->outward[1]),point->outward[2]);
        point->valid=point->degree==2 && nl>0 && ol>0 && isfinite(nl) && isfinite(ol);
        if(point->valid)for(int d=0;d<3;d++){point->normal[d]/=nl;point->outward[d]/=ol;}
    }
    return front;
}

static int asr_fit(Arena_T arena,const AsmChart *a,const AsmChart *b,
                   const AsmCorr *samples,const double *weights,size_t count,int parity,AsmRelation *relation)
{
    /* The same two source-face operators used by the repair field. An
     * orthonormal-frame gap is not itself an intrinsic UV displacement;
     * carry it through its original face before fitting either direction.
     * Midpoint fitting can cancel their disagreement and understate RMS. */
    if(count>SIZE_MAX/(4*sizeof(double)))return 0;
    double *qa=ARENA_ALLOC(arena,4*count*sizeof(double)),*qb=ARENA_ALLOC(arena,4*count*sizeof(double));
    double *mass=ARENA_ALLOC(arena,2*count*sizeof(double));
    AsmChart local_a=*a,local_b=*b;
    local_a.flags&=~ASM_CHART_MIRROR;local_b.flags&=~ASM_CHART_MIRROR;
    local_a.placed_uv=local_b.placed_uv=NULL;
    local_a.pose_x=local_a.pose_y=local_a.pose_theta=0;
    local_b.pose_x=local_b.pose_y=local_b.pose_theta=0;
    for(size_t k=0;k<count;k++){
        const AsmCorr *c=samples+k;
        if(c->valid!=3 || !c->trim_a.valid || !c->trim_b.valid || !(weights[k]>0) || !isfinite(weights[k]))return 0;
        if(!act_trim_point(&local_a,c->va,&c->trim_a,c->gap_a,1,qa+4*k) ||
           !act_trim_point(&local_b,c->vb,&c->trim_b,c->gap_b,-1,qb+4*k+2))return 0;
        for(int d=0;d<2;d++){
            qa[4*k+2+d]=a->uv[2*(size_t)c->va+d];qb[4*k+d]=b->uv[2*(size_t)c->vb+d];
        }
        if(parity){qb[4*k]=-qb[4*k];qb[4*k+2]=-qb[4*k+2];}
        mass[2*k]=mass[2*k+1]=.5*weights[k];
    }
    double theta,tx,ty,rms;
    if(AsmRelate_rigid_fit(qa,qb,mass,2*count,&theta,&tx,&ty,&rms) ||
       !isfinite(theta) || !isfinite(tx) || !isfinite(ty) || !isfinite(rms))return 0;
    relation->theta=theta;relation->tx=tx;relation->ty=ty;relation->rms=rms;
    return 1;
}

/* The pass reads chart pairs whose cubes share an origin or a face
 * (asr_domain), and a pair's outcome depends only on its two charts. So the
 * candidates come from a per-origin chart index instead of all n^2 pairs, the
 * existing-relation test is a binary search instead of a scan of every
 * relation, the independent per-pair work runs in parallel, and outcomes are
 * published serially in (ai, bi) order: the relations, correspondences and
 * their order are exactly those of the original nested loop. */
typedef struct AsrCand {
    int32_t a, b, axis, status;   /* status: 0 skip, 1 relation, -1 error */
    double plane;
    AsmRelation relation;         /* valid when status == 1; corr_first set on publish */
    AsmCorr *samples;             /* status == 1: count samples in a per-thread keep arena */
    size_t count;
} AsrCand;
typedef struct AsrOrigin { long z, y, x; int32_t cube; } AsrOrigin;

static int asr_origin_order(const void *p, const void *q)
{
    const AsrOrigin *x = p, *y = q;
    if (x->z != y->z) return x->z < y->z ? -1 : 1;
    if (x->y != y->y) return x->y < y->y ? -1 : 1;
    if (x->x != y->x) return x->x < y->x ? -1 : 1;
    return x->cube < y->cube ? -1 : x->cube > y->cube;
}
static int asr_int64_order(const void *p, const void *q)
{
    int64_t x = *(const int64_t *)p, y = *(const int64_t *)q;
    return x < y ? -1 : x > y;
}
static int asr_int32_order(const void *p, const void *q)
{
    int32_t x = *(const int32_t *)p, y = *(const int32_t *)q;
    return x < y ? -1 : x > y;
}
static int64_t asr_pair_key(int32_t a, int32_t b)
{ return a < b ? ((int64_t)a << 32) | (uint32_t)b : ((int64_t)b << 32) | (uint32_t)a; }

/* [first, last) of the origins equal to (z, y, x) in the sorted table */
static size_t asr_origin_range(const AsrOrigin *o, size_t n, long z, long y, long x, size_t *last)
{
    size_t lo = 0, hi = n;
    AsrOrigin key = { z, y, x, INT32_MIN };
    while (lo < hi) { size_t mid = lo + (hi - lo) / 2; if (asr_origin_order(o + mid, &key) < 0) lo = mid + 1; else hi = mid; }
    size_t end = lo;
    while (end < n && o[end].z == z && o[end].y == y && o[end].x == x) end++;
    *last = end; return lo;
}

/* One pair's reconstruction: the body of the original loop. Returns -1 for
 * the cases that stopped the original pass, 0 for the ones it skipped, and 1
 * with the relation and a copy of its samples in `keep`. */
static int asr_pair(Arena_T scratch, Arena_T keep, const AsmChart *a, const AsmChart *b, AsrCand *out)
{
    int axis = out->axis; double plane = out->plane;
    Arena_Mark mark = Arena_save(scratch); int status = 0;
    uint8_t *va=ARENA_CALLOC(scratch,a->nv,1),*vb=ARENA_CALLOC(scratch,b->nv,1);
    ActFrame *fa=act_vertex_frames(scratch,a,va),*fb=act_vertex_frames(scratch,b,vb);
    AsmTrim *ta=act_trim_bases(scratch,a,fa,va),*tb=act_trim_bases(scratch,b,fb,vb);
    if(!ta || !tb){status=-1;goto done;}
    AsrFront *front_a=NULL,*front_b=NULL;
    if(axis<0){
        front_a=asr_front_geometry(scratch,a);front_b=asr_front_geometry(scratch,b);
        if(!front_a || !front_b){status=-1;goto done;}
    }
    int32_t *ia=ARENA_ALLOC(scratch,a->nv*sizeof(int32_t)),*ib=ARENA_ALLOC(scratch,b->nv*sizeof(int32_t));
    float *xa=ARENA_ALLOC(scratch,3*a->nv*sizeof(float)),*xb=ARENA_ALLOC(scratch,3*b->nv*sizeof(float));
    size_t na=0,nb=0;
    for(size_t v=0;v<a->nv;v++)if(a->boundary[v] && va[v] && (axis<0?front_a[v].valid:fabs(a->xyz[3*v+axis]-plane)<=2)){ia[na]=(int32_t)v;memcpy(xa+3*na++,a->xyz+3*v,3*sizeof(float));}
    for(size_t v=0;v<b->nv;v++)if(b->boundary[v] && vb[v] && (axis<0?front_b[v].valid:fabs(b->xyz[3*v+axis]-plane)<=2)){ib[nb]=(int32_t)v;memcpy(xb+3*nb++,b->xyz+3*v,3*sizeof(float));}
    if(na<6 || nb<6)goto done;
    KDTree_T ka=KDTree_new(scratch,xa,na),kb=KDTree_new(scratch,xb,nb);
    AsmCorr *samples=ARENA_CALLOC(scratch,na,sizeof(AsmCorr));size_t count=0;int parity=-1,inconsistent=0;
    double sumgap=0,sumnormal=0;
    for(size_t k=0;k<na;k++){
        float dist;size_t j=KDTree_nearest(kb,xa+3*k,&dist);
        if(dist>16 || KDTree_nearest(ka,xb+3*j,NULL)!=k)continue;
        double delta[3],an[3],bn[3],dot=0,ag=0,bg=0,frame_dot=0;
        asr_normal(fa+ia[k],an);asr_normal(fb+ib[j],bn);
        for(int d=0;d<3;d++)frame_dot+=an[d]*bn[d];
        if(axis<0)for(int d=0;d<3;d++){an[d]=front_a[ia[k]].normal[d];bn[d]=front_b[ib[j]].normal[d];}
        for(int d=0;d<3;d++){delta[d]=(double)xb[3*j+d]-xa[3*k+d];dot+=an[d]*bn[d];ag+=delta[d]*an[d];bg+=delta[d]*bn[d];}
        if(fabs(dot)<.75 || fmax(fabs(ag),fabs(bg))>1)continue;
        if(axis<0){
            double opposed=0,toward_a=0,toward_b=0;
            for(int d=0;d<3;d++){
                double oa=front_a[ia[k]].outward[d],ob=front_b[ib[j]].outward[d];
                opposed-=oa*ob;toward_a+=oa*delta[d];toward_b-=ob*delta[d];
            }
            if(!(opposed>0) || toward_a<0 || toward_b<0)continue;
        }
        int p=frame_dot<0;if(parity<0)parity=p;else if(parity!=p)inconsistent=1;
        AsmCorr c={0};c.va=ia[k];c.vb=ib[j];c.run=-1;c.trim_a=ta[c.va];c.trim_b=tb[c.vb];
        if(!c.trim_a.valid || !c.trim_b.valid || !act_frame_gap(fa+c.va,delta,c.gap_a) || !act_frame_gap(fb+c.vb,delta,c.gap_b)){inconsistent=1;continue;}
        c.valid=3;samples[count++]=c;sumgap+=sqrt(dist);sumnormal+=fmax(fabs(ag),fabs(bg));
    }
    AsmBoundaryStats ordered={0};
    if(count<6 || inconsistent || AsmBoundary_order(scratch,a,b,samples,count,&ordered) || !ordered.supported_runs)goto done;
    double *weights=ARENA_ALLOC(scratch,count*sizeof(double));
    if(AsmBoundary_weights(scratch,a,b,samples,count,weights)){status=-1;goto done;}
    double mass=0;for(size_t k=0;k<count;k++)mass+=weights[k];
    if(!(mass>0) || !isfinite(mass)){status=-1;goto done;}
    AsmRelation relation={0};relation.a=out->a;relation.b=out->b;
    relation.flags=ASM_REL_WEAK|(parity?ASM_REL_PARITY:0);relation.continuity=ASM_CONT_SOURCE|ASM_CONT_RECONSTRUCTED;
    if(!asr_fit(scratch,a,b,samples,weights,count,parity,&relation)){status=-1;goto done;}
    relation.corr_count=(int32_t)(count<=(size_t)INT32_MAX?count:0);relation.n_corr=count;
    relation.seam_len=mass;relation.mean_gap=sumgap/count;relation.normal_gap=sumnormal/count;relation.robust_w=1;
    relation.w_xy=mass;relation.w_theta=mass;
    out->relation=relation;out->count=count;
    out->samples=ARENA_ALLOC(keep,count*sizeof(AsmCorr));memcpy(out->samples,samples,count*sizeof(AsmCorr));
    status=1;
done:
    Arena_restore(scratch,mark);
    return status;
}

int AsmContinuity_reconstruct(AsmRun *run)
{
    if(!run || !run->arena)return -1;
    if(!run->pile || !run->n_cubes)return 0;
    size_t nc=run->n_charts;
    Arena_T scratch=Arena_new();int rc=0;size_t added=0,internal_added=0;
    /* the existing (min, max) relation pairs at entry; this pass visits each
     * unordered pair once, so its own additions can never be looked up */
    size_t nkeys=run->n_rels;
    int64_t *keys=ARENA_ALLOC(scratch,(nkeys?nkeys:1)*sizeof(int64_t));
    for(size_t r=0;r<nkeys;r++)keys[r]=asr_pair_key(run->rels[r].a,run->rels[r].b);
    qsort(keys,nkeys,sizeof *keys,asr_int64_order);
    /* charts of each identified cube origin, in ascending chart order */
    AsrOrigin *origin=ARENA_ALLOC(scratch,run->n_cubes*sizeof *origin);size_t n_origin=0;
    for(size_t c=0;c<run->n_cubes;c++)if(run->pile[c].has_id){
        origin[n_origin].z=(long)run->pile[c].oz;origin[n_origin].y=(long)run->pile[c].oy;origin[n_origin].x=(long)run->pile[c].ox;
        origin[n_origin].cube=(int32_t)c;n_origin++;
    }
    qsort(origin,n_origin,sizeof *origin,asr_origin_order);
    size_t *first=ARENA_CALLOC(scratch,run->n_cubes+1,sizeof(size_t));
    for(size_t i=0;i<nc;i++){int32_t q=run->charts[i].cube;if(q>=0 && (size_t)q<run->n_cubes)first[q+1]++;}
    for(size_t q=0;q<run->n_cubes;q++)first[q+1]+=first[q];
    int32_t *members=ARENA_ALLOC(scratch,(nc?nc:1)*sizeof(int32_t));size_t *fill=ARENA_CALLOC(scratch,run->n_cubes,sizeof(size_t));
    for(size_t i=0;i<nc;i++){int32_t q=run->charts[i].cube;if(q>=0 && (size_t)q<run->n_cubes)members[first[q]+fill[q]++]=(int32_t)i;}
    long step=(long)llround(run->cube_size);int neighbours=run->cube_size>0 && fabs((double)step-run->cube_size)<=1e-8;
    /* the candidates, filtered exactly as the original loop filtered them
     * before its per-pair work, in (ai, bi) order */
    size_t cap=1024,ncand=0;AsrCand *cand=malloc(cap*sizeof *cand);
    int32_t *partner=ARENA_ALLOC(scratch,(nc?nc:1)*sizeof(int32_t));
    if(!cand){Arena_dispose(&scratch);return -1;}
    for(size_t ai=0;ai<nc;ai++){
        const AsmChart *a=run->charts+ai;
        if(!AsmChart_in_layout(a) || !a->boundary || a->cube<0 || (size_t)a->cube>=run->n_cubes || !run->pile[a->cube].has_id)continue;
        const MeshPileEntry *pa=run->pile+a->cube;size_t np=0;
        long oz=(long)pa->oz,oy=(long)pa->oy,ox=(long)pa->ox;
        long look[7][3]={{oz,oy,ox},{oz-step,oy,ox},{oz+step,oy,ox},{oz,oy-step,ox},{oz,oy+step,ox},{oz,oy,ox-step},{oz,oy,ox+step}};
        for(int l=0;l<(neighbours?7:1);l++){
            size_t last,o=asr_origin_range(origin,n_origin,look[l][0],look[l][1],look[l][2],&last);
            for(;o<last;o++){
                int32_t q=origin[o].cube;
                for(size_t m=first[q];m<first[q+1];m++)if((size_t)members[m]>ai)partner[np++]=members[m];
            }
        }
        qsort(partner,np,sizeof *partner,asr_int32_order);
        for(size_t k=0;k<np;k++){
            if(k && partner[k]==partner[k-1])continue;
            size_t bi=(size_t)partner[k];const AsmChart *b=run->charts+bi;int axis;double plane;
            if(!AsmChart_in_layout(b) || !b->boundary || !asr_domain(run,a,b,&axis,&plane))continue;
            int64_t key=asr_pair_key((int32_t)ai,(int32_t)bi);
            if(bsearch(&key,keys,nkeys,sizeof *keys,asr_int64_order))continue; /* a rejected measured proposal is still evidence */
            int separated=0;for(int d=0;d<3;d++)if(a->bbox_lo[d]>b->bbox_hi[d]+4 || b->bbox_lo[d]>a->bbox_hi[d]+4)separated=1;
            if(separated)continue;
            if(ncand==cap){AsrCand *grown=realloc(cand,2*cap*sizeof *cand);if(!grown){free(cand);Arena_dispose(&scratch);return -1;}cand=grown;cap*=2;}
            AsrCand *c=cand+ncand++;memset(c,0,sizeof *c);
            c->a=(int32_t)ai;c->b=(int32_t)bi;c->axis=axis;c->plane=plane;
        }
    }
    int threads=1;
#ifdef _OPENMP
    threads=omp_get_max_threads();if(threads<1)threads=1;
#endif
    Arena_T *work=calloc((size_t)threads,sizeof *work),*keep=calloc((size_t)threads,sizeof *keep);
    if(!work || !keep){free(work);free(keep);free(cand);Arena_dispose(&scratch);return -1;}
    for(int t=0;t<threads;t++){work[t]=Arena_new();keep[t]=Arena_new();}
    int n_int=ncand<=(size_t)INT_MAX?(int)ncand:-1,k=0;
    if(n_int<0)rc=-1;
    else{
#pragma omp parallel for schedule(dynamic,8) num_threads(threads) if(n_int>=64)
        for(k=0;k<n_int;k++){
            int t=0;
#ifdef _OPENMP
            t=omp_get_thread_num();
#endif
            cand[k].status=asr_pair(work[t],keep[t],run->charts+cand[k].a,run->charts+cand[k].b,cand+k);
        }
    }
    /* publish in the original order; stop where the original stopped */
    for(size_t i=0;i<ncand && !rc;i++){
        AsrCand *c=cand+i;
        if(c->status<0){rc=-1;break;}
        if(c->status==0)continue;
        if(c->count>INT32_MAX || run->n_corr>INT32_MAX-c->count){rc=-1;break;}
        c->relation.corr_first=(int32_t)run->n_corr;c->relation.corr_count=(int32_t)c->count;
        for(size_t s=0;s<c->count;s++)AsmRun_push_corr(run,c->samples+s);
        AsmRun_push_rel(run,&c->relation);added++;internal_added+=c->axis<0;
    }
    fprintf(stderr,"[assemble source] reconstructed %zu omitted chart-pair boundary hypotheses (%zu internal facing fronts); source identity remains unreviewed\n",added,internal_added);
    fprintf(stderr,"[assemble source] %zu same-cube or face-adjacent candidate pairs searched on %d threads\n",ncand,threads);
    for(int t=0;t<threads;t++){Arena_dispose(&work[t]);Arena_dispose(&keep[t]);}
    free(work);free(keep);free(cand);Arena_dispose(&scratch);return rc;
}
