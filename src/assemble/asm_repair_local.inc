/* Direct scalar form of unroller_v2/contact_projection.py: unique shared
 * vertices, twelve source-triangle axes, exact one-dimensional LS projection.
 * Included by asm_repair.c so the implementation shares its solver lifetime. */
static void apr_rotation(const double j[4], double r[4])
{
    double c=j[0]+j[3],s=j[2]-j[1],length=hypot(c,s);
    if(length>0){c/=length;s/=length;}else{c=1;s=0;}
    r[0]=c;r[1]=-s;r[2]=s;r[3]=c;
}

static double apr_project_contact(const int32_t fa[3],const int32_t fb[3],
                                  const double *uv,AprContactCopy *copy,double delta[12])
{
    int32_t ids[6];for(int k=0;k<3;k++){ids[k]=fa[k];ids[k+3]=fb[k];}
    for(int i=1;i<6;i++){int32_t id=ids[i];int j=i;while(j && ids[j-1]>id){ids[j]=ids[j-1];j--;}ids[j]=id;}
    copy->count=0;
    for(int k=0;k<6;k++)if(!k || ids[k]!=ids[k-1])copy->ids[copy->count++]=ids[k];
    int ina[6]={0},inb[6]={0};double points[12],center[2]={0};
    for(int i=0;i<copy->count;i++){
        for(int k=0;k<3;k++){ina[i]|=copy->ids[i]==fa[k];inb[i]|=copy->ids[i]==fb[k];}
        for(int d=0;d<2;d++)center[d]+=(uv[2*(size_t)copy->ids[i]+d]-uv[2*(size_t)copy->ids[0]+d])/copy->count;
    }
    for(int i=0;i<copy->count;i++)for(int d=0;d<2;d++)points[2*i+d]=uv[2*(size_t)copy->ids[i]+d]-uv[2*(size_t)copy->ids[0]+d]-center[d];
    double best=INFINITY;memset(delta,0,12*sizeof(double));
    for(int axis=0;axis<12;axis++){
        const int32_t *face=axis%6<3?fa:fb;int k=axis%3;double sign=axis<6?1:-1;
        double nx=sign*(uv[2*(size_t)face[(k+1)%3]+1]-uv[2*(size_t)face[k]+1]);
        double ny=sign*(uv[2*(size_t)face[k]]-uv[2*(size_t)face[(k+1)%3]]);
        double length=hypot(nx,ny);if(!(length>0))continue;nx/=length;ny/=length;
        double scalar[6],sum=0,square=0;int order[6],count=0;
        for(int i=0;i<copy->count;i++){
            scalar[i]=nx*points[2*i]+ny*points[2*i+1];order[i]=i;
            if(inb[i]){count++;sum+=scalar[i];square+=scalar[i]*scalar[i];}
        }
        for(int i=1;i<copy->count;i++){int at=order[i],j=i;while(j && scalar[order[j-1]]>scalar[at]){order[j]=order[j-1];j--;}order[j]=at;}
        double height=0,cost=INFINITY;
        for(int interval=0;interval<=copy->count;interval++){
            double lo=interval?scalar[order[interval-1]]:-INFINITY;
            double hi=interval<copy->count?scalar[order[interval]]:INFINITY;
            double h=fmax(lo,fmin(hi,count?sum/count:0));
            double value=count*h*h-2*sum*h+square;
            if(isfinite(h) && value<cost){height=h;cost=value;}
            if(interval<copy->count){int at=order[interval],role=ina[at]-inb[at];count+=role;sum+=role*scalar[at];square+=role*scalar[at]*scalar[at];}
        }
        double shifts[6],actual=0;
        for(int i=0;i<copy->count;i++){
            double corrected=ina[i]&&inb[i]?height:ina[i]?fmax(scalar[i],height):fmin(scalar[i],height);
            shifts[i]=corrected-scalar[i];actual+=shifts[i]*shifts[i];
        }
        if(actual<best){best=actual;for(int i=0;i<copy->count;i++){delta[2*i]=shifts[i]*nx;delta[2*i+1]=shifts[i]*ny;}}
    }
    return best;
}

static void apr_contact_copy(void *context,size_t a,size_t b,double area)
{
    AprSolver *s=context;AsmField *f=s->field;AprContactCopy copy={0};double delta[12];
    apr_preserve_self_contact(context,a,b,area);
    if(s->opts.contact_copy_active &&
       !s->opts.contact_copy_active(f,a,b,s->opts.contact_merit_context))return;
    double cost=apr_project_contact(f->faces+3*a,f->faces+3*b,s->evaluated_uv,&copy,delta);
    copy.weight=s->opts.contact_weight*fmin(f->metric[a].area,f->metric[b].area)/
                (3*s->opts.contact_length*s->opts.contact_length);
    s->separation_energy+=copy.weight*cost;
    if(!s->evaluated_gradient)return;
    if(s->local_count==s->local_capacity){
        size_t next=s->local_capacity?2*s->local_capacity:256;
        if(next<s->local_capacity || next>SIZE_MAX/sizeof(AprContactCopy)){s->allocation_failed=1;return;}
        AprContactCopy *p=realloc(s->local_contacts,next*sizeof *p);
        if(!p){s->allocation_failed=1;return;}s->local_contacts=p;s->local_capacity=next;
    }
    s->local_contacts[s->local_count++]=copy;
    for(int k=0;k<copy.count;k++){
        size_t v=(size_t)copy.ids[k];
        for(int d=0;d<2;d++)s->gradient[2*v+d]-=2*copy.weight*delta[2*k+d];
        s->diagonal[v]+=2*copy.weight*(1.-1./copy.count);
    }
}
