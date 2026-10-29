/* New material has no inherited contact allowance. Existing neighbors may
 * retain contacts only within the same old chart pair, without more area.
 * Triangle subdivision cannot turn a shrinking overlap into a regression.
 * Self contacts retain the stronger original
 * face-pair check. The caller still checks the complete external field. */
typedef struct ApaContactGroup {
    uint64_t key;
    size_t before_pairs,after_pairs;
    double before_area,after_area,self_regression_area;
    int used;
} ApaContactGroup;

/* A refusal by a handful of new-material pairs is traced pair by pair
 * (residual_contact): the last sliver a solve could not clear. */
enum {APA_RESIDUAL_PAIRS=8};

typedef struct ApaContactCertificate {
    Arena_T arena;
    const AsmField *field;
    const uint8_t *core;
    const double *reference;
    ApaContactGroup *groups;
    size_t capacity,count,new_pairs,self_regressions;
    double new_area;
    int before,okay;
    size_t residual[2*APA_RESIDUAL_PAIRS];
    double residual_area[APA_RESIDUAL_PAIRS];
} ApaContactCertificate;

static double apa_face_gap(const AsmField *f,size_t a,size_t b);

/* 1 when every corner of triangle `in` lies inside (or on) triangle `out`. */
static int apa_triangle_inside(const double out[6],const double in[6])
{
    double det=(out[2]-out[0])*(out[5]-out[1])-(out[4]-out[0])*(out[3]-out[1]);
    if(!(fabs(det)>0))return 0;
    for(int k=0;k<3;k++){
        double x=in[2*k]-out[0],y=in[2*k+1]-out[1];
        double s=(x*(out[5]-out[1])-(out[4]-out[0])*y)/det,t=((out[2]-out[0])*y-x*(out[3]-out[1]))/det;
        if(s<-1e-12 || t<-1e-12 || s+t>1+1e-12)return 0;
    }
    return 1;
}

static void apa_residual_trace(const AsmField *f,const uint8_t *moving,const ApaContactCertificate *s,FILE *trace)
{
    size_t n=s->new_pairs<APA_RESIDUAL_PAIRS?s->new_pairs:APA_RESIDUAL_PAIRS;
    for(size_t k=0;k<n;k++){
        size_t face[2]={s->residual[2*k],s->residual[2*k+1]};double tri[2][6],g[2][6],norm[2]={0,0};int edge[2][3];
        double area=AsmContacts_pair_faces(f->contacts,f->uv,face[0],face[1],g[0],g[1]);
        for(int side=0;side<2;side++){
            size_t c=(size_t)f->face_chart[face[side]];const AsmChart *ch=f->run->charts+c;
            for(int j=0;j<3;j++){
                size_t v=(size_t)f->faces[3*face[side]+(size_t)j];
                tri[side][2*j]=f->uv[2*v];tri[side][2*j+1]=f->uv[2*v+1];
                size_t local=v-f->vertex_offset[c];
                edge[side][j]=ch->boundary && local<ch->nv?ch->boundary[local]:-1;
                norm[side]+=g[side][2*j]*g[side][2*j]+g[side][2*j+1]*g[side][2*j+1];
            }
        }
        int inside=apa_triangle_inside(tri[0],tri[1])?1:apa_triangle_inside(tri[1],tri[0])?2:0;
        fprintf(trace,"{\"phase\":\"residual_contact\",\"pair\":%zu,\"area\":%.17g,\"recomputed_area\":%.17g,\"gap3d\":%.17g,"
            "\"contained\":\"%s\"",k,s->residual_area[k],area,apa_face_gap(f,face[0],face[1]),
            inside==1?"b_in_a":inside==2?"a_in_b":"partial");
        for(int side=0;side<2;side++){
            size_t c=(size_t)f->face_chart[face[side]];
            fprintf(trace,",\"%c\":{\"chart\":%zu,\"face\":%zu,\"core\":%d,\"moving\":%d,\"boundary\":[%d,%d,%d],"
                "\"gradient_norm\":%.17g,\"uv\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]}",side?'b':'a',c,
                face[side]-f->face_offset[c],s->core[c]!=0,moving[c]!=0,edge[side][0],edge[side][1],edge[side][2],sqrt(norm[side]),
                tri[side][0],tri[side][1],tri[side][2],tri[side][3],tri[side][4],tri[side][5]);
        }
        fputs("}\n",trace);
    }
}

static size_t apa_contact_slot(uint64_t key,size_t capacity)
{
    key^=key>>33;key*=UINT64_C(0xff51afd7ed558ccd);key^=key>>33;
    return (size_t)key&(capacity-1);
}

static ApaContactGroup *apa_contact_group(ApaContactCertificate *s,uint64_t key)
{
    if(!s->okay)return NULL;
    if(!s->capacity || s->count>=s->capacity/2){
        if(s->capacity>SIZE_MAX/(2*sizeof *s->groups)){s->okay=0;return NULL;}
        size_t capacity=s->capacity?2*s->capacity:64;
        ApaContactGroup *groups=ARENA_CALLOC(s->arena,capacity,sizeof *groups);
        for(size_t i=0;i<s->capacity;i++)if(s->groups[i].used){
            size_t slot=apa_contact_slot(s->groups[i].key,capacity);
            while(groups[slot].used)slot=(slot+1)&(capacity-1);
            groups[slot]=s->groups[i];
        }
        s->groups=groups;s->capacity=capacity;
    }
    size_t slot=apa_contact_slot(key,s->capacity);
    while(s->groups[slot].used && s->groups[slot].key!=key)slot=(slot+1)&(s->capacity-1);
    ApaContactGroup *group=s->groups+slot;
    if(!group->used){group->used=1;group->key=key;s->count++;}
    return group;
}

static void apa_contact_visit(void *context,size_t a,size_t b,double area)
{
    ApaContactCertificate *s=context;const AsmField *f=s->field;
    uint32_t ca=(uint32_t)f->face_chart[a],cb=(uint32_t)f->face_chart[b];
    if(s->core[ca] || s->core[cb]){
        if(!s->before){
            if(s->new_pairs<APA_RESIDUAL_PAIRS){
                s->residual[2*s->new_pairs]=a;s->residual[2*s->new_pairs+1]=b;s->residual_area[s->new_pairs]=area;
            }
            s->new_pairs++;s->new_area+=area;
        }
        return;
    }
    uint64_t key=((uint64_t)(ca<cb?ca:cb)<<32)|(uint64_t)(ca<cb?cb:ca);
    ApaContactGroup *group=apa_contact_group(s,key);if(!group)return;
    if(s->before){group->before_pairs++;group->before_area+=area;return;}
    group->after_pairs++;group->after_area+=area;
    if(ca==cb){
        double old=AsmContacts_pair_faces(f->contacts,s->reference,a,b,NULL,NULL);
        if(!(old>0) || area>old+64*DBL_EPSILON*fmax(1,old)){
            s->self_regressions++;group->self_regression_area+=area;
        }
    }
}

/* Scheduling observes the same obligations as admission. An unchanged old
 * overlap is not unfinished new-material work. New chart-pair and self-contact
 * regressions still contribute even when aggregate old area decreases. */
static double apa_blocking_contact_area(const ApaContactCertificate *s)
{
    double area=s->new_area;
    for(size_t i=0;i<s->capacity;i++)if(s->groups[i].used){
        const ApaContactGroup *g=s->groups+i;
        double excess=g->after_area>g->before_area+1e-9*fmax(1,g->before_area)?g->after_area-g->before_area:0;
        area+=!g->before_pairs && g->after_pairs?g->after_area:fmax(excess,g->self_regression_area);
    }
    return area;
}

static const char *apa_contact_certificate(AsmField *f,const uint8_t *moving,const uint8_t *core,
                                          const double *reference,const uint8_t *active,FILE *trace,int *complete)
{
    if(complete)*complete=0;
    int neighbors=0;for(size_t c=0;c<f->run->n_charts;c++)neighbors|=moving[c] && !core[c];
    AsmContactStats contacts={0};AsmContacts_active_faces(f->contacts,active);
    ApaContactCertificate certificate={0};certificate.arena=f->arena;certificate.field=f;
    certificate.core=core;certificate.reference=reference;certificate.okay=1;
    if(neighbors){
        certificate.before=1;
        int rc=AsmContacts_measure(f->contacts,reference,NULL,apa_contact_visit,&certificate,&contacts);
        if(rc || !contacts.complete || !certificate.okay)
            return rc==-2?"contact_walk_over_budget":"contact_query_incomplete";
        certificate.before=0;
    }
    /* Without moving neighbors every pair is new material; the visit then
     * only records pairs for the residual trace. */
    int record=neighbors || trace;
    int rc=AsmContacts_measure(f->contacts,f->uv,NULL,record?apa_contact_visit:NULL,
                              record?&certificate:NULL,&contacts);
    if(!neighbors){certificate.new_pairs=contacts.pairs;certificate.new_area=contacts.area;certificate.self_regressions=0;}
    size_t regressions=0;
    if(neighbors)for(size_t i=0;i<certificate.capacity;i++)if(certificate.groups[i].used){
        const ApaContactGroup *g=certificate.groups+i;
        regressions+=(!g->before_pairs && g->after_pairs) ||
            g->after_area>g->before_area+1e-9*fmax(1,g->before_area);
    }
    int finished=!rc && contacts.complete && certificate.okay;
    if(complete)*complete=finished;
    if(trace)fprintf(trace,"{\"phase\":\"final_contacts\",\"pairs\":%zu,\"area\":%.17g,\"complete\":%s,"
        "\"new_material_pairs\":%zu,\"new_material_area\":%.17g,\"neighbor_pairs\":%zu,"
        "\"neighbor_chart_pair_regressions\":%zu,\"self_contact_regressions\":%zu,"
        "\"neighbor_contact_policy\":\"nonincreasing_chart_pair_area_v2\"}\n",
        contacts.pairs,contacts.area,finished?"true":"false",certificate.new_pairs,certificate.new_area,
        contacts.pairs-certificate.new_pairs,regressions,certificate.self_regressions);
    if(!finished)return rc==-2?"contact_walk_over_budget":"contact_query_incomplete";
    if(certificate.new_pairs && trace && certificate.new_pairs<=APA_RESIDUAL_PAIRS)
        apa_residual_trace(f,moving,&certificate,trace);
    if(certificate.new_pairs)return "new_material_contact_failed";
    return regressions || certificate.self_regressions?"incumbent_contact_regression":NULL;
}

/* Reuse the immutable incumbent allowance while measuring intermediate
 * candidates. These are local scores, never replacement certificates; the
 * final admission still checks the entire outside sheet and original faces. */
typedef struct ApaSolverContacts {
    ApaContactCertificate certificate;
    ApaContactCertificate raw;
    FILE *trace;
    size_t evaluations;
    int ready,raw_ready;
} ApaSolverContacts;

/* A continuous proposal objective: new material pays for its full overlap;
 * inherited chart pairs pay only for area above their original raw overlap.
 * Full boundary contributions survive, including tolerated material contact.
 * This surrogate never replaces the material-level admission certificate. */
static double apa_raw_contact_collect(void *context,size_t a,size_t b,double area)
{
    ApaContactCertificate *s=context;const AsmField *f=s->field;
    uint32_t ca=(uint32_t)f->face_chart[a],cb=(uint32_t)f->face_chart[b];
    if(s->core[ca] || s->core[cb]){if(!s->before)s->new_area+=area;return 0;}
    uint64_t key=((uint64_t)(ca<cb?ca:cb)<<32)|(uint64_t)(ca<cb?cb:ca);
    ApaContactGroup *g=apa_contact_group(s,key);if(!g)return NAN;
    if(s->before)g->before_area+=area;else g->after_area+=area;
    return 0;
}

static int apa_solver_force_active(const AsmField *f,size_t a,size_t b,void *context)
{
    ApaSolverContacts *state=context;ApaContactCertificate *s=&state->raw;
    uint32_t ca=(uint32_t)f->face_chart[a],cb=(uint32_t)f->face_chart[b];
    if(s->core[ca] || s->core[cb])return 1;
    uint64_t key=((uint64_t)(ca<cb?ca:cb)<<32)|(uint64_t)(ca<cb?cb:ca);
    ApaContactGroup *g=apa_contact_group(s,key);
    return g && g->after_area>g->before_area;
}

static double apa_raw_contact_weight(void *context,size_t a,size_t b,double area)
{
    (void)area;ApaSolverContacts *state=context;
    return (double)apa_solver_force_active(state->raw.field,a,b,context);
}

static int apa_solver_copy_active(const AsmField *f,size_t a,size_t b,void *context)
{
    ApaSolverContacts *state=context;ApaContactCertificate *s=&state->raw;
    uint32_t ca=(uint32_t)f->face_chart[a],cb=(uint32_t)f->face_chart[b];
    if(s->core[ca] || s->core[cb])return 1;
    uint64_t key=((uint64_t)(ca<cb?ca:cb)<<32)|(uint64_t)(ca<cb?cb:ca);
    ApaContactGroup *g=apa_contact_group(s,key);
    /* Raw excess supplies the inherited-pair force. Switching on its whole
     * old projection cost at an infinitesimal increase would create a jump
     * in energy at the inherited allowance instead of a continuous penalty. */
    return g && !(g->before_area>0);
}

static int apa_solver_contact_objective(const AsmField *f,const double *uv,double *gradient,
    AsmContactVisit visit,void *visit_context,AsmContactStats *contacts,double *energy,void *context)
{
    ApaSolverContacts *state=context;ApaContactCertificate *raw=&state->raw;
    AsmContactStats collected={0};double ignored=0;
    if(state->raw_ready<0)return -1;
    if(!state->raw_ready){
        raw->arena=state->certificate.arena;raw->field=f;raw->core=state->certificate.core;
        raw->reference=state->certificate.reference;raw->okay=1;raw->before=1;
        int rc=AsmContacts_weighted_objective(f->contacts,raw->reference,NULL,NULL,NULL,
            apa_raw_contact_collect,raw,&collected,&ignored);
        raw->before=0;
        if(rc || !collected.complete || !raw->okay){state->raw_ready=-1;return rc?rc:-1;}
        state->raw_ready=1;
    }
    raw->new_area=0;
    for(size_t k=0;k<raw->capacity;k++)raw->groups[k].after_area=0;
    int rc=AsmContacts_weighted_objective(f->contacts,uv,NULL,NULL,NULL,
        apa_raw_contact_collect,raw,&collected,&ignored);
    if(rc || !collected.complete || !raw->okay)return rc?rc:-1;
    /* Select whole chart pairs after their totals are known. Selecting only
     * certified triangle intersections would erase boundary forces again. */
    double value=raw->new_area;
    for(size_t k=0;k<raw->capacity;k++)if(raw->groups[k].used)
        value+=fmax(0,raw->groups[k].after_area-raw->groups[k].before_area);
    rc=AsmContacts_weighted_objective(f->contacts,uv,gradient,visit,visit_context,
        apa_raw_contact_weight,state,contacts,&ignored);
    if(rc || !contacts->complete || !raw->okay)return rc?rc:-1;
    *energy=value;return 0;
}

static int apa_solver_contact_merit(const AsmField *field,const double *uv,double *area,void *context)
{
    ApaSolverContacts *state=context;ApaContactCertificate *c=&state->certificate;
    AsmContactStats contacts={0};
    if(state->ready<0)return -1;
    if(!state->ready){
        c->before=1;
        int rc=AsmContacts_measure(field->contacts,c->reference,NULL,apa_contact_visit,c,&contacts);
        c->before=0;
        if(rc || !contacts.complete || !c->okay){state->ready=-1;return -1;}
        state->ready=1;
    }
    c->new_pairs=c->self_regressions=0;c->new_area=0;
    for(size_t k=0;k<c->capacity;k++){
        c->groups[k].after_pairs=0;c->groups[k].after_area=c->groups[k].self_regression_area=0;
    }
    int rc=AsmContacts_measure(field->contacts,uv,NULL,apa_contact_visit,c,&contacts);
    if(rc || !contacts.complete || !c->okay)return -1;
    *area=apa_blocking_contact_area(c);
    if(state->trace)fprintf(state->trace,"{\"phase\":\"blocking_contact_progress\",\"evaluation\":%zu,"
        "\"new_material_area\":%.17g,\"inherited_regression_area\":%.17g,\"blocking_area\":%.17g,"
        "\"total_area\":%.17g,\"new_material_pairs\":%zu,\"self_regressions\":%zu}\n",
        ++state->evaluations,c->new_area,*area-c->new_area,*area,contacts.area,c->new_pairs,c->self_regressions);
    return 0;
}
