typedef struct SRGaugeEvent { double angle, db, dc; } SRGaugeEvent;

static int sr_gauge_compare(const void *left,const void *right)
{
    double a=((const SRGaugeEvent *)left)->angle,b=((const SRGaugeEvent *)right)->angle;
    return a<b?-1:a>b;
}

static void sr_gauge_interval(double b,double c,double lo,double hi,
                               double *best,double *angle)
{
    double candidates[3]={lo,hi,atan2(c,b)};
    for (int k=0;k<3;k++) if (candidates[k]>=lo && candidates[k]<=hi) {
        double value=b*cos(candidates[k])+c*sin(candidates[k]);
        if (value>*best) { *best=value;*angle=candidates[k]; }
    }
}

/* For one common original-frame turn, each face's cross energy is
 * weight*(constant-|b*cos(2*t)+c*sin(2*t)|). Its sign-change angles partition
 * one quarter turn into intervals with an exact sinusoidal optimum. */
static int sr_direction_gauge_proposal(Arena_T arena,const SRSystem *system,
                                       const double *current,const double *original,
                                       double *proposal,double *angle)
{
    const double pi=3.14159265358979323846;
    Arena_Mark mark=Arena_save(arena);
    SRGaugeEvent *events=ARENA_ALLOC(arena,system->direction_faces*sizeof *events);
    int8_t *sign=ARENA_ALLOC(arena,system->ncomponents*sizeof *sign);
    for (size_t c=0;c<system->ncomponents;c++) sign[c]=1;
    for (size_t f=0;f<system->nf;f++) {
        const int32_t *face=system->faces+3*f;
        sign[system->vertex_component[face[0]]]=sr_orient2(original+2*(size_t)face[0],
            original+2*(size_t)face[1],original+2*(size_t)face[2])<0?-1:1;
    }
    size_t count=0;double sum_b=0,sum_c=0;
    for (size_t f=0;f<system->nf;f++) if (system->direction[f].weight>0) {
        const SRDirection *d=system->direction+f;const int32_t *face=system->faces+3*f;
        double h[2],v[2];sr_direction_vectors(d,face,current,h,v);
        double b=.5*(h[1]*h[1]+v[0]*v[0]-h[0]*h[0]-v[1]*v[1]);
        double c=(h[0]*h[1]-v[0]*v[1])*sign[system->vertex_component[face[0]]];
        double s=b>0?1:b<0?-1:c>0?1:c<0?-1:0;
        sum_b+=d->weight*s*b;sum_c+=d->weight*s*c;
        double crossing=fmod(atan2(-b,c)+pi,pi);
        if (crossing>0 && crossing<pi && s) {
            events[count].angle=crossing;
            events[count].db=-2*d->weight*s*b;events[count].dc=-2*d->weight*s*c;count++;
        }
    }
    if (!isfinite(sum_b) || !isfinite(sum_c)) { Arena_restore(arena,mark);return -1; }
    qsort(events,count,sizeof *events,sr_gauge_compare);
    double previous=0,best=sum_b,best_angle=0;
    for (size_t e=0;e<count;) {
        double next=events[e].angle;
        sr_gauge_interval(sum_b,sum_c,previous,next,&best,&best_angle);
        do { sum_b+=events[e].db;sum_c+=events[e].dc;e++; } while (e<count && events[e].angle==next);
        previous=next;
    }
    sr_gauge_interval(sum_b,sum_c,previous,pi,&best,&best_angle);
    *angle=.5*best_angle;if (*angle>.25*pi) *angle-=.5*pi;
    double center[2]={0,0};
    for (size_t v=0;v<system->nv;v++) {
        center[0]+=current[2*v]/(double)system->nv;
        center[1]+=current[2*v+1]*sign[system->vertex_component[v]]/(double)system->nv;
    }
    double co=cos(*angle),si=sin(*angle);
    for (size_t v=0;v<system->nv;v++) {
        double s=sign[system->vertex_component[v]];
        double x=current[2*v]-center[0],y=s*current[2*v+1]-center[1];
        proposal[2*v]=center[0]+co*x-si*y;
        proposal[2*v+1]=s*(center[1]+si*x+co*y);
    }
    Arena_restore(arena,mark);return 0;
}
