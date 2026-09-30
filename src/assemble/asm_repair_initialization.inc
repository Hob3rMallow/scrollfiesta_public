/* assembly_initialization.translated_hull_seed: discrete whole-chart motion
 * from original triangle contact regions, followed by complete geometry
 * checks. These are temporary seeds, never final seam certificates. */
typedef struct AhiContact {int32_t a,b;size_t fa,fb;} AhiContact;
typedef struct AhiPoint {double x,y;} AhiPoint;
typedef struct AhiPlane {int32_t a,b;size_t na,nb;int32_t *vertices;double normal[2],required;} AhiPlane;
typedef struct AhiRegions {AprSolver *s;int *column;AhiContact *pairs;size_t n,capacity;int failed;} AhiRegions;

static void ahi_contact(void *context,size_t fa,size_t fb,double area)
{
    AhiRegions *r=context;const AsmField *f=r->s->field;int32_t a=f->face_chart[fa],b=f->face_chart[fb];(void)area;
    if(r->failed || a==b || (r->column[a]<0 && r->column[b]<0))return;
    if(a>b){int32_t c=a;a=b;b=c;size_t t=fa;fa=fb;fb=t;}
    if(r->n==r->capacity){
        if(r->capacity>SIZE_MAX/(2*sizeof(AhiContact))){r->failed=1;return;}
        size_t capacity=r->capacity?2*r->capacity:128;AhiContact *next=realloc(r->pairs,capacity*sizeof *next);
        if(!next){r->failed=1;return;}r->pairs=next;r->capacity=capacity;
    }
    r->pairs[r->n++]=(AhiContact){a,b,fa,fb};
}
static int ahi_pair_order(const void *a,const void *b)
{
    const AhiContact *x=a,*y=b;if(x->a!=y->a)return(x->a>y->a)-(x->a<y->a);
    return(x->b>y->b)-(x->b<y->b);
}
static int ahi_point_order(const void *a,const void *b)
{
    const AhiPoint *x=a,*y=b;if(x->x!=y->x)return(x->x>y->x)-(x->x<y->x);
    return(x->y>y->y)-(x->y<y->y);
}
static double ahi_cross(AhiPoint a,AhiPoint b,AhiPoint c)
{return(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);}
static AhiPoint *ahi_hull(Arena_T arena,const double *uv,const int32_t *ids,size_t n,size_t *size)
{
    *size=0;if(n<3 || n>SIZE_MAX/(2*sizeof(AhiPoint)))return NULL;
    AhiPoint *points=ARENA_ALLOC(arena,n*sizeof *points),*hull=ARENA_ALLOC(arena,2*n*sizeof *hull);
    for(size_t i=0;i<n;i++)points[i]=(AhiPoint){uv[2*(size_t)ids[i]],uv[2*(size_t)ids[i]+1]};
    qsort(points,n,sizeof *points,ahi_point_order);size_t unique=0;
    for(size_t i=0;i<n;i++)if(!unique || points[i].x!=points[unique-1].x || points[i].y!=points[unique-1].y)points[unique++]=points[i];
    if(unique<3)return NULL;
    size_t k=0;
    for(size_t i=0;i<unique;i++){while(k>=2 && ahi_cross(hull[k-2],hull[k-1],points[i])<=0)k--;hull[k++]=points[i];}
    size_t split=k+1;
    for(size_t i=unique-1;i-- >0;){while(k>=split && ahi_cross(hull[k-2],hull[k-1],points[i])<=0)k--;hull[k++]=points[i];}
    if(k<4)return NULL;*size=k-1;return hull;
}
static int ahi_separation(const AhiPoint *a,size_t na,const AhiPoint *b,size_t nb,double normal[2],double *gap)
{
    double direction[2]={0};for(size_t i=0;i<na;i++){direction[0]-=a[i].x/na;direction[1]-=a[i].y/na;}
    for(size_t i=0;i<nb;i++){direction[0]+=b[i].x/nb;direction[1]+=b[i].y/nb;}
    *gap=INFINITY;
    for(int side=0;side<2;side++){
        const AhiPoint *hull=side?b:a;size_t n=side?nb:na;
        for(size_t i=0;i<n;i++){
            size_t next=(i+1)%n;double nx=hull[i].y-hull[next].y,ny=hull[next].x-hull[i].x,length=hypot(nx,ny);
            if(!(length>1e-12))continue;nx/=length;ny/=length;
            if(nx*direction[0]+ny*direction[1]<0){nx=-nx;ny=-ny;}
            double high=-INFINITY,low=INFINITY;
            for(size_t k=0;k<na;k++)high=fmax(high,nx*a[k].x+ny*a[k].y);
            for(size_t k=0;k<nb;k++)low=fmin(low,nx*b[k].x+ny*b[k].y);
            double value=high-low+.5;if(value<*gap){*gap=value;normal[0]=nx;normal[1]=ny;}
        }
    }
    return isfinite(*gap);
}
static size_t ahi_unique(int32_t *ids,size_t n)
{
    qsort(ids,n,sizeof *ids,ars_int_order);size_t count=0;
    for(size_t i=0;i<n;i++)if(!count || ids[i]!=ids[count-1])ids[count++]=ids[i];return count;
}
static int ahi_add_regions(ArsProblem *p,AhiRegions *r,const double *reference,AhiPlane **planes,size_t *count,size_t *capacity)
{
    const AsmField *f=p->s->field;AsmContactStats contacts={0};r->n=0;
    if(AsmContacts_touching(f->contacts,f->uv,ahi_contact,r,&contacts) || r->failed)return -1;
    if(r->n>1)qsort(r->pairs,r->n,sizeof *r->pairs,ahi_pair_order);int fresh=0;
    for(size_t first=0;first<r->n;){
        size_t end=first+1;while(end<r->n && r->pairs[end].a==r->pairs[first].a && r->pairs[end].b==r->pairs[first].b)end++;
        size_t n=3*(end-first);if(n>SIZE_MAX/(2*sizeof(int32_t)))return -1;
        Arena_T temporary=Arena_new();int32_t *va=ARENA_ALLOC(temporary,2*n*sizeof(int32_t)),*vb=va+n;
        for(size_t k=first;k<end;k++){memcpy(va+3*(k-first),f->faces+3*r->pairs[k].fa,3*sizeof(int32_t));memcpy(vb+3*(k-first),f->faces+3*r->pairs[k].fb,3*sizeof(int32_t));}
        size_t na=ahi_unique(va,n),nb=ahi_unique(vb,n);int32_t a=r->pairs[first].a,b=r->pairs[first].b;int found=0;
        for(size_t k=0;k<*count;k++){
            const AhiPlane *plane=*planes+k;
            if(plane->a==a && plane->b==b && plane->na==na && plane->nb==nb &&
                !memcmp(plane->vertices,va,na*sizeof(int32_t)) && !memcmp(plane->vertices+na,vb,nb*sizeof(int32_t))){found=1;break;}
        }
        if(!found){
            size_t ha,hb;AhiPoint *one=ahi_hull(temporary,f->uv,va,na,&ha),*two=ahi_hull(temporary,f->uv,vb,nb,&hb);
            double normal[2],gap;if(!one || !two || !ahi_separation(one,ha,two,hb,normal,&gap)){Arena_dispose(&temporary);return -1;}
            if(gap>1e-7){
                if(*count==*capacity){
                    if(*capacity>SIZE_MAX/(2*sizeof(AhiPlane))){Arena_dispose(&temporary);return -1;}
                    size_t next_capacity=*capacity?2* *capacity:64;AhiPlane *next=realloc(*planes,next_capacity*sizeof *next);
                    if(!next){Arena_dispose(&temporary);return -1;}*planes=next;*capacity=next_capacity;
                }
                double high=-INFINITY,low=INFINITY;
                for(size_t k=0;k<na;k++)high=fmax(high,normal[0]*reference[2*(size_t)va[k]]+normal[1]*reference[2*(size_t)va[k]+1]);
                for(size_t k=0;k<nb;k++)low=fmin(low,normal[0]*reference[2*(size_t)vb[k]]+normal[1]*reference[2*(size_t)vb[k]+1]);
                int32_t *vertices=ARENA_ALLOC(p->arena,(na+nb)*sizeof(int32_t));memcpy(vertices,va,na*sizeof(int32_t));memcpy(vertices+na,vb,nb*sizeof(int32_t));
                (*planes)[(*count)++]=(AhiPlane){a,b,na,nb,vertices,{normal[0],normal[1]},high-low+.5};fresh++;
            }
        }
        Arena_dispose(&temporary);first=end;
    }
    return fresh;
}
/* 1 separated, 0 refused or bounded, -1 malformed/allocation failure. A caller
 * must keep its incumbent separately: a translation seed may open seams. */
static int apr_initialize_separated(AprSolver *s,int passes)
{
    if(apr_expired(&s->opts))return 0;
    AsmField *f=s->field;ArsProblem p={0};p.s=s;p.arena=Arena_new();AhiRegions regions={0};regions.s=s;
    AhiPlane *planes=NULL;size_t count=0,capacity=0;int rc=0;size_t nc=f->run->n_charts,bytes=2*f->nv*sizeof(double);
    int *column=ARENA_ALLOC(p.arena,nc*sizeof(int));regions.column=column;
    for(size_t c=0;c<nc;c++){column[c]=-1;if(f->vertex_offset[c]!=SIZE_MAX && !s->rigid_fixed[c]){
        if(p.n>INT_MAX-2){rc=-1;goto done;}column[c]=p.n;p.n+=2;
    }}
    if(!p.n)goto done;
    double *reference=ARENA_ALLOC(p.arena,bytes);memcpy(reference,f->uv,bytes);
    double *target=ARENA_CALLOC(p.arena,(size_t)p.n,sizeof(double));p.answer=ARENA_CALLOC(p.arena,(size_t)p.n,sizeof(double));
    double *area=ARENA_CALLOC(p.arena,nc,sizeof(double));for(size_t t=0;t<f->nf;t++)area[f->face_chart[t]]+=f->metric[t].area;
    for(size_t c=0;c<nc;c++)if(column[c]>=0)for(int d=0;d<2;d++)alq_entry(&p.h,column[c]+d,column[c]+d,4*area[c]/4096);
    for(size_t si=0;si<f->n_seams;si++)if(ars_incident(s,f->seams+si)){
        const AsmFieldSeam *seam=f->seams+si;
        for(size_t k=0;k<seam->count;k++)for(int side=0;side<2;side++){
            const AsmWitness *w=f->witness+seam->first+k;if(!w->valid || !(w->weight>0)){rc=-1;goto done;}
            double residual[2],coefficient[2]={0};AsmField_residual(w,side,reference,residual);
            int charts[]={seam->a,seam->b};
            for(int q=0;q<4;q++)coefficient[f->vertex_chart[w->vertex[side][q]]==seam->b]+=w->coefficient[side][q];
            for(int a=0;a<2;a++)if(column[charts[a]]>=0)for(int d=0;d<2;d++){
                int row=column[charts[a]]+d;target[row]-=32*w->weight*coefficient[a]*residual[d];
                for(int b=0;b<2;b++)if(column[charts[b]]>=0 && row>=column[charts[b]]+d)
                    alq_entry(&p.h,row,column[charts[b]]+d,32*w->weight*coefficient[a]*coefficient[b]);
            }
        }
    }
    if(p.h.failed){rc=-1;goto done;}
    for(int iteration=0;iteration<passes && !apr_expired(&s->opts);iteration++){
        int fresh=ahi_add_regions(&p,&regions,reference,&planes,&count,&capacity);
        if(fresh<0){rc=-1;goto done;}
        if(!fresh){
            rc=1;for(size_t k=0;k<count;k++){
                const AhiPlane *plane=planes+k;double delta[2]={0};
                for(int d=0;d<2;d++)delta[d]=(column[plane->b]>=0?p.answer[column[plane->b]+d]:0)-(column[plane->a]>=0?p.answer[column[plane->a]+d]:0);
                if(plane->normal[0]*delta[0]+plane->normal[1]*delta[1]<plane->required-1e-6)rc=0;
            }
            break;
        }
        p.m=0;p.a.n=0;p.a.failed=0;
        for(size_t k=0;k<count;k++){
            const AhiPlane *plane=planes+k;int ids[4],n=0;double values[4];
            for(int side=0;side<2;side++){
                int c=side?plane->b:plane->a;if(column[c]<0)continue;
                for(int d=0;d<2;d++){ids[n]=column[c]+d;values[n++]=(side?-1:1)*plane->normal[d];}
            }
            if(!ars_plane(&p,n,ids,values,-plane->required-1e-5)){rc=-1;goto done;}
        }
        s->stats->initialization_trials++;s->stats->initialization_constraints+=(size_t)fresh;
        if(!ConvexQp_translation_seed(p.n,p.h.n,p.h.row,p.h.col,p.h.value,target,p.m,p.a.n,p.a.row,p.a.col,p.a.value,p.upper,p.answer)){
            fprintf(stderr,"[assemble initialization] iteration %d: %d new regions, %zu constraints; QP refused\n",iteration+1,fresh,count);break;
        }
        s->stats->initialization_steps++;double maximum=0;
        for(size_t c=0;c<nc;c++)if(column[c]>=0){
            maximum=fmax(maximum,hypot(p.answer[column[c]],p.answer[column[c]+1]));
            for(size_t k=0;k<f->run->charts[c].nv;k++){size_t v=f->vertex_offset[c]+k;
                for(int d=0;d<2;d++)f->uv[2*v+d]=reference[2*v+d]+p.answer[column[c]+d];}
        }
        /* These are discrete chart placements, not a small deformation of
         * the old hierarchy. Refit alone can mix distant chart fragments in
         * every leaf and turn the next complete query into an all-pairs scan. */
        if(AsmContacts_rebuild(f->contacts,f->uv)){rc=-1;goto done;}
        fprintf(stderr,"[assemble initialization] iteration %d: %d new regions, %zu constraints, maximum translation %.9g\n",iteration+1,fresh,count,maximum);
    }
    if(rc==1){
        AsmContactStats contacts={0};if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&contacts)){rc=-1;goto done;}
        if(contacts.pairs)rc=0;
        for(size_t at=0;at<s->n_active && rc==1;at++){
            size_t t=s->active_faces[at];double j[4],lo,hi,det;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
            if(!UvGuard_positive(f->uv,f->faces+3*t) || !(lo>0) || !isfinite(hi))rc=0;
        }
    }
done:
    s->stats->initialization_separated+=rc==1;free(regions.pairs);free(planes);alq_free(&p.h);alq_free(&p.a);free(p.upper);Arena_dispose(&p.arena);
    return rc;
}

/* The same maximum-support forest used by Python's rigid_initialize and C
 * admission. Every chart containing a fixed/gauge vertex remains a root;
 * cycles stay required. Failure rolls back every temporary chart transform. */
static int apr_restart_forest(AprSolver *s)
{
    if(!s->opts.iterations || !s->n_active || apr_expired(&s->opts))return 0;
    AsmField *f=s->field;Arena_T arena=Arena_new();size_t nc=f->run->n_charts,bytes=2*f->nv*sizeof(double);
    uint8_t *core=ARENA_CALLOC(arena,nc,1),*known=ARENA_CALLOC(arena,nc,1);
    double *entry=ARENA_ALLOC(arena,bytes);memcpy(entry,f->uv,bytes);
    for(size_t at=0;at<s->n_active;at++)core[f->face_chart[s->active_faces[at]]]=1;
    for(size_t c=0;c<nc;c++)known[c]=(uint8_t)(f->vertex_offset[c]!=SIZE_MAX && s->rigid_fixed[c]);
    s->stats->forest_trials++;int ok=apa_initialize(f,core,known);
    for(size_t v=0;v<f->nv && ok;v++){
        if(!isfinite(f->uv[2*v]) || !isfinite(f->uv[2*v+1]))ok=0;
        if(s->fixed[v] && memcmp(f->uv+2*v,entry+2*v,2*sizeof(double)))ok=0;
    }
    if(!ok)memcpy(f->uv,entry,bytes);
    int rc=ok;
    if(AsmContacts_rebuild(f->contacts,f->uv)){
        memcpy(f->uv,entry,bytes);(void)AsmContacts_rebuild(f->contacts,f->uv);rc=-1;
    }
    if(rc>0)s->stats->forest_initialized++;
    fprintf(stderr,"[assemble initialization] maximum-support forest: %s; all source cycles retained\n",rc<0?"contact index failed":rc?"initialized":"refused");
    Arena_dispose(&arena);return rc;
}

/* A failed warm ordering says nothing about physical feasibility. Retry a
 * source-forest placement from the immutable entry, then separate its actual
 * original contact regions. The caller owns the guarded incumbent archive. */
static int apr_initialize_retry(AprSolver *s,const double *entry,int force_forest)
{
    if(apr_expired(&s->opts))return 0;
    if(!force_forest){int rc=apr_initialize_separated(s,24);if(rc!=0)return rc;}
    memcpy(s->field->uv,entry,2*s->field->nv*sizeof(double));
    int rc=apr_restart_forest(s);if(rc>0)rc=apr_initialize_separated(s,24);
    if(rc || !s->initial_uv || !memcmp(s->initial_uv,entry,2*s->field->nv*sizeof(double)))return rc;
    /* Earlier local/global deformation can obstruct a forest that separates
     * the original shapes. Retain that alternative throughout this solve.
     * The caller's current incumbent and acceptance gates remain unchanged. */
    for(size_t v=0;v<s->field->nv;v++)if(s->fixed[v] &&
        memcmp(s->initial_uv+2*v,entry+2*v,2*sizeof(double)))return 0;
    memcpy(s->field->uv,s->initial_uv,2*s->field->nv*sizeof(double));
    s->stats->forest_entry_restarts++;
    fprintf(stderr,"[assemble initialization] retrying original solver-entry shapes; current incumbent retained\n");
    rc=apr_restart_forest(s);return rc>0?apr_initialize_separated(s,24):rc;
}
