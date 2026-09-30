#include "asm_reading_order.h"
#include "asm_store.h"
#include "asm_continuity.h"
#include "asm_field.h"
#include "../whole/atlas_winding_sync.h"
#include "../common/union_find.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ARO_TAU 6.283185307179586476925286766559

typedef struct AroChart {
    size_t first,count,piece,group,island;
    double theta,axial_lo,axial_hi,radius,pitch,qlo,qhi;
    double angular_gradient[2];
    int32_t lift;
    int valid,phase_conflict,relative;
} AroChart;
typedef struct AroEdge { int32_t a,b,delta; } AroEdge;
typedef struct AroModel {
    AroChart *chart;
    size_t charts,samples,pieces,source_edges,phase_conflicts;
    size_t relative_charts,within_charts,unresolved_charts,local_pitch_charts;
    double relative_area,within_area,unresolved_area,pitch;
    AtlasWindingSyncStats sync;
    int radial_solved;
    int sense;
    size_t sense_pieces;
    AtlasWindingSyncRelation *relations;
    size_t nrelations,radial_cycle_conflicts;
} AroModel;

static double aro_wrap(double x) { return x-ARO_TAU*floor(x/ARO_TAU+0.5); }
static double aro_area(const AsmChart *ch,size_t f)
{
    const int32_t *t=ch->faces+3*f;double u[3],v[3],n[3];
    for(int d=0;d<3;d++){u[d]=(double)ch->xyz[3*(size_t)t[1]+d]-ch->xyz[3*(size_t)t[0]+d];v[d]=(double)ch->xyz[3*(size_t)t[2]+d]-ch->xyz[3*(size_t)t[0]+d];}
    n[0]=u[1]*v[2]-u[2]*v[1];n[1]=u[2]*v[0]-u[0]*v[2];n[2]=u[0]*v[1]-u[1]*v[0];
    return .5*sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
}
static size_t aro_count(const AsmChart *ch)
{
    if(!AsmChart_registered(ch) || !ch->nf || !(ch->area3d>0))return 0;
    double n=ceil(ch->area3d/ASM_READING_SAMPLE_AREA);
    size_t count=(size_t)fmin(ASM_READING_MAX_SAMPLES,fmax(ASM_READING_MIN_SAMPLES,n));
    return count<ch->nf?count:ch->nf;
}
static double aro_vertex_phase(const AsmChart *ch,int32_t v,const AsmAxis *axis,double center)
{
    double xyz[3],p[3];for(int d=0;d<3;d++)xyz[d]=ch->xyz[3*(size_t)v+d];
    AsmAxis_frame(axis,xyz,p);return center+aro_wrap(atan2(p[2],p[1])-center);
}
static int aro_source_edge(const AsmRun *run,const AsmAxis *axis,const AroChart *chart,
                           const AsmSeamMeasure *certificate,size_t i,AroEdge *out)
{
    const AsmRelation *r=run->rels+i;
    if(!chart[r->a].valid || !chart[r->b].valid ||
       (r->flags&(ASM_REL_DROPPED|ASM_REL_CONTACT|ASM_REL_CROSSWRAP)) ||
       !(r->continuity&(ASM_CONT_REQUIRED|ASM_CONT_VERIFIED|ASM_CONT_CUT)))return 0;
    if(!certificate[i].pass)return 0;
    int count=0,delta=0;
    for(int k=0;k<3;k++){
        int j=k*(r->corr_count-1)/2;if(j<0)continue;
        const AsmCorr *w=run->corr+r->corr_first+j;
        if(w->valid!=3 || w->va<0 || w->vb<0 || (size_t)w->va>=run->charts[r->a].nv || (size_t)w->vb>=run->charts[r->b].nv)continue;
        double a=aro_vertex_phase(run->charts+r->a,w->va,axis,chart[r->a].theta);
        double b=aro_vertex_phase(run->charts+r->b,w->vb,axis,chart[r->b].theta);
        double gap=(a-b)/ARO_TAU;int current=(int)lround(gap);
        if(fabs(gap-current)>ASM_READING_PHASE_RESIDUAL || (count && current!=delta))return 0;
        delta=current;count++;
    }
    if(!count)return 0;out->a=r->a;out->b=r->b;out->delta=delta;return 1;
}

/* Integer lifting is exact graph traversal, followed by ALL-edge cycle
 * checks. It supplies relative phase only, not an absolute winding prior. */
static size_t aro_lift(Arena_T arena,AroChart *chart,size_t nc,const AroEdge *edges,size_t ne,size_t *conflicts)
{
    size_t *offset=ARENA_CALLOC(arena,nc+1,sizeof *offset),*cursor=ARENA_ALLOC(arena,nc*sizeof *cursor);
    for(size_t e=0;e<ne;e++){offset[edges[e].a+1]++;offset[edges[e].b+1]++;}
    for(size_t c=1;c<=nc;c++)offset[c]+=offset[c-1];memcpy(cursor,offset,nc*sizeof *cursor);
    size_t *adj=ARENA_ALLOC(arena,(2*ne+1)*sizeof *adj),*queue=ARENA_ALLOC(arena,(nc+1)*sizeof *queue);
    for(size_t e=0;e<ne;e++){adj[cursor[edges[e].a]++]=e;adj[cursor[edges[e].b]++]=e;}
    for(size_t c=0;c<nc;c++)chart[c].piece=SIZE_MAX;
    size_t pieces=0;
    for(size_t c=0;c<nc;c++)if(chart[c].valid && chart[c].piece==SIZE_MAX){
        size_t head=0,tail=0;chart[c].piece=pieces++;chart[c].lift=0;queue[tail++]=c;
        while(head<tail){size_t a=queue[head++];
            for(size_t k=offset[a];k<offset[a+1];k++){
                const AroEdge *e=edges+adj[k];size_t b=(size_t)(e->a==(int32_t)a?e->b:e->a);
                if(chart[b].piece!=SIZE_MAX)continue;
                chart[b].piece=chart[a].piece;chart[b].lift=chart[a].lift+(e->a==(int32_t)a?e->delta:-e->delta);queue[tail++]=b;
            }
        }
    }
    uint8_t *bad=ARENA_CALLOC(arena,pieces+1,1);*conflicts=0;
    for(size_t e=0;e<ne;e++)if(chart[edges[e].b].lift-chart[edges[e].a].lift!=edges[e].delta){bad[chart[edges[e].a].piece]=1;(*conflicts)++;}
    for(size_t c=0;c<nc;c++)if(chart[c].valid)chart[c].phase_conflict=bad[chart[c].piece];
    return pieces;
}

/* Determine spiral handedness from long source-connected pieces. Regress
 * radius on relative phase while removing linear axial taper; no absolute
 * radius-to-wind rounding is used. Short pieces cannot determine the sign. */
static int aro_sense(Arena_T arena,const AroChart *chart,size_t nc,size_t np,
                     const float *phi,const double *radius,const double *axial,double pitch,size_t *used)
{
    double *mom=ARENA_CALLOC(arena,11*(np+1),sizeof *mom);
    for(size_t p=0;p<np;p++){mom[11*p+9]=DBL_MAX;mom[11*p+10]=-DBL_MAX;}
    for(size_t c=0;c<nc;c++)if(chart[c].valid && !chart[c].phase_conflict){
        const AroChart *ch=chart+c;double *s=mom+11*ch->piece;
        for(size_t f=ch->first;f<ch->first+ch->count;f++){
            double q=((double)phi[3*f]+phi[3*f+1]+phi[3*f+2])/(3*ARO_TAU),z=axial[f],r=radius[f];
            s[0]++;s[1]+=q;s[2]+=z;s[3]+=r;s[4]+=q*q;s[5]+=z*z;s[6]+=q*z;s[7]+=q*r;s[8]+=z*r;
            s[9]=fmin(s[9],q);s[10]=fmax(s[10],q);
        }
    }
    double vote=0;*used=0;
    for(size_t p=0;p<np;p++){
        double *s=mom+11*p,n=s[0];if(n<3 || s[10]-s[9]<ASM_READING_SENSE_MIN_TURNS)continue;
        double qq=s[4]-s[1]*s[1]/n,zz=s[5]-s[2]*s[2]/n,qz=s[6]-s[1]*s[2]/n;
        double qr=s[7]-s[1]*s[3]/n,zr=s[8]-s[2]*s[3]/n,den=qq*zz-qz*qz;
        if(!(den>1e-8*fmax(1,qq*zz)))continue;
        double slope=(qr*zz-zr*qz)/den;
        if(fabs(slope)<ASM_READING_SENSE_MIN_PITCH*pitch || fabs(slope)>ASM_READING_SENSE_MAX_PITCH*pitch)continue;
        vote+=slope>0?n:-n;(*used)++;
    }
    return vote<0?-1:1;
}

static int aro_model_mode(Arena_T arena,const AsmRun *run,const AsmAxis *axis,AroModel *out,int radial)
{
    memset(out,0,sizeof *out);size_t nc=run->n_charts,capacity=0;
    if(nc>(size_t)INT32_MAX)return -1;
    AroChart *chart=ARENA_CALLOC(arena,nc+1,sizeof *chart);out->chart=chart;out->charts=nc;
    out->pitch=run->layer_d_global;
    for(size_t c=0;c<nc;c++)capacity+=aro_count(run->charts+c);
    if(capacity>(size_t)INT32_MAX/3)return -1;
    float *vertices=ARENA_ALLOC(arena,(9*capacity+1)*sizeof *vertices),*phi=ARENA_ALLOC(arena,(3*capacity+1)*sizeof *phi);
    int32_t *faces=ARENA_ALLOC(arena,(3*capacity+1)*sizeof *faces),*face_chart=ARENA_ALLOC(arena,(capacity+1)*sizeof *face_chart);
    int32_t *islands=ARENA_ALLOC(arena,(capacity+1)*sizeof *islands);
    double *radius=ARENA_ALLOC(arena,(capacity+1)*sizeof *radius),*axial=ARENA_ALLOC(arena,(capacity+1)*sizeof *axial),*pitch=ARENA_ALLOC(arena,(capacity+1)*sizeof *pitch);
    size_t ns=0;
    for(size_t c=0;c<nc;c++){
        const AsmChart *ch=run->charts+c;size_t wanted=aro_count(ch);AroChart *q=chart+c;
        q->piece=q->group=q->island=SIZE_MAX;if(!wanted)continue;
        double center[3];AsmAxis_frame(axis,ch->centroid,center);q->theta=atan2(center[2],center[1]);
        q->pitch=run->chart_layer_d?run->chart_layer_d[c]:0;
        if(!isfinite(q->pitch) || q->pitch<0)q->pitch=0;
        out->local_pitch_charts+=q->pitch>0;
        q->first=ns;q->axial_lo=q->qlo=DBL_MAX;q->axial_hi=q->qhi=-DBL_MAX;
        double cumulative=0,next=.5*ch->area3d/(double)wanted,min_radius=DBL_MAX;
        for(size_t f=0;f<ch->nf && q->count<wanted;f++){
            cumulative+=aro_area(ch,f);if(cumulative<next && f+1<ch->nf)continue;
            double rsum=0,ssum=0;
            for(int k=0;k<3;k++){
                size_t v=(size_t)ch->faces[3*f+k];double xyz[3],frame[3];
                for(int d=0;d<3;d++)vertices[9*ns+3*k+d]=(float)(xyz[d]=ch->xyz[3*v+d]);
                AsmAxis_frame(axis,xyz,frame);double r=hypot(frame[1],frame[2]);min_radius=fmin(min_radius,r);
                double theta=q->theta+aro_wrap(atan2(frame[2],frame[1])-q->theta);
                phi[3*ns+k]=(float)theta;faces[3*ns+k]=(int32_t)(3*ns+k);rsum+=r;ssum+=frame[0];
                q->qlo=fmin(q->qlo,theta/ARO_TAU);q->qhi=fmax(q->qhi,theta/ARO_TAU);
                q->axial_lo=fmin(q->axial_lo,frame[0]);q->axial_hi=fmax(q->axial_hi,frame[0]);
            }
            radius[ns]=rsum/3;axial[ns]=ssum/3;pitch[ns]=q->pitch;face_chart[ns]=(int32_t)c;q->radius+=radius[ns];
            /* A sampled face represents equal source area. Its local phase
             * gradient is branch independent and chooses the ribbon's U
             * direction without an elastic deformation or reflection. */
            if(!radial){
                double uv[3][2];for(int k=0;k<3;k++)AsmChart_point(ch,(size_t)ch->faces[3*f+k],uv[k]);
                double u=uv[1][0]-uv[0][0],v=uv[1][1]-uv[0][1],a=uv[2][0]-uv[0][0],b=uv[2][1]-uv[0][1];
                double det=u*b-v*a,dq=aro_wrap(phi[3*ns+1]-phi[3*ns]),eq=aro_wrap(phi[3*ns+2]-phi[3*ns]);
                if(fabs(det)>1e-20){double weight=ch->area3d/(double)wanted;
                    q->angular_gradient[0]+=weight*(dq*b-eq*v)/det;
                    q->angular_gradient[1]+=weight*(eq*u-dq*a)/det;}
            }
            ns++;q->count++;next=((double)q->count+.5)*ch->area3d/(double)wanted;
        }
        if(!q->count || min_radius<ASM_READING_CORE_RADIUS || q->qhi-q->qlo>.5){ns=q->first;q->count=0;continue;}
        q->radius/=(double)q->count;q->valid=1;
        if(c && !(c%5000))fprintf(stderr,"[reading order] sampled %zu/%zu charts, %zu original triangles\n",c,nc,ns);
    }
    out->samples=ns;
    AroEdge *edge=ARENA_ALLOC(arena,(run->n_rels+1)*sizeof *edge);size_t ne=0;
    AsmSeamMeasure *certificate=ARENA_ALLOC(arena,(run->n_rels+1)*sizeof *certificate);
    if(AsmField_measure_source_seams(run,certificate))return -1;
    for(size_t e=0;e<run->n_rels;e++)if(aro_source_edge(run,axis,chart,certificate,e,edge+ne))ne++;
    out->source_edges=ne;size_t np=aro_lift(arena,chart,nc,edge,ne,&out->phase_conflicts);out->pieces=np;
    for(size_t c=0;c<nc;c++)if(chart[c].valid){
        AroChart *q=chart+c;q->qlo+=q->lift;q->qhi+=q->lift;
        for(size_t f=q->first;f<q->first+q->count;f++){
            islands[f]=(int32_t)q->piece;face_chart[f]=(int32_t)q->piece;
            for(int k=0;k<3;k++)phi[3*f+k]+=(float)(ARO_TAU*q->lift);
        }
    }
    out->sense=aro_sense(arena,chart,nc,np,phi,radius,axial,out->pitch,&out->sense_pieces);
    if(out->sense<0)for(size_t c=0;c<nc;c++)if(chart[c].valid){double low=chart[c].qlo;chart[c].qlo=-chart[c].qhi;chart[c].qhi=-low;}
    fprintf(stderr,"[reading order] %zu source pieces, %zu cycle conflicts; spiral sense %d from %zu long pieces\n",np,out->phase_conflicts,out->sense,out->sense_pieces);
    int32_t *fi=NULL,*correction=NULL;size_t ni=0,nr=0;double *prior=NULL,*mad=NULL;AtlasWindingSyncRelation *relations=NULL;
    AtlasWindingSyncProblem problem={0};problem.vertices=vertices;problem.nvertices=3*ns;problem.faces=faces;problem.nfaces=ns;
    problem.phi=phi;problem.face_radius=radius;problem.face_axial=axial;problem.face_chart=face_chart;problem.ncharts=np;
    problem.pitch=out->pitch;problem.spiral_b=out->pitch/ARO_TAU;problem.sense=out->sense;
    problem.axial_bin_spacing=ASM_READING_AXIAL_BIN;problem.phase_bins=ASM_READING_PHASE_BINS;
    problem.source_face_island=islands;problem.nsource_islands=np;problem.face_pitch=pitch;
    out->radial_solved=radial && ns && out->pitch>0 && !AtlasWindingSync_solve(arena,&problem,&fi,&ni,&correction,&prior,&mad,&relations,&nr,&out->sync);
    out->relations=relations;out->nrelations=nr;
    /* The legacy radius prior is useful for its original atlas, but is not
     * an identity observation. Relative reading order uses exact integer
     * differences and tests every eligible cycle, with a free group gauge. */
    AroChart *relative=ARENA_CALLOC(arena,np+1,sizeof *relative);
    if(out->radial_solved){
        AroEdge *radial_edges=ARENA_ALLOC(arena,(nr+1)*sizeof *radial_edges);size_t count=0;
        for(size_t p=0;p<ni;p++)relative[p].valid=1;
        for(size_t r=0;r<nr;r++)if(relations[r].eligible)
            radial_edges[count++]=(AroEdge){relations[r].island0,relations[r].island1,relations[r].target_turn_correction};
        aro_lift(arena,relative,ni,radial_edges,count,&out->radial_cycle_conflicts);
        out->sync.satisfied_relations=0;
        for(size_t r=0;r<nr;r++){
            AtlasWindingSyncRelation *e=relations+r;
            e->solved_turn_correction=relative[e->island1].lift-relative[e->island0].lift;
            e->final_residual=e->solved_turn_correction-e->target_turn_correction;
            out->sync.satisfied_relations+=e->eligible && !e->final_residual;
        }
    }
    size_t *piece_island=ARENA_ALLOC(arena,(np+1)*sizeof *piece_island);
    for(size_t p=0;p<np;p++)piece_island[p]=p;
    if(out->radial_solved)for(size_t f=0;f<ns;f++)piece_island[islands[f]]=(size_t)fi[f];
    UnionFind components=UF_new(arena,(int32_t)(np?np:1));
    uint8_t *bad=ARENA_CALLOC(arena,np+1,1),*supported=ARENA_CALLOC(arena,np+1,1);
    if(out->radial_solved){
        for(size_t r=0;r<nr;r++)if(relations[r].eligible)uf_union(&components,relations[r].island0,relations[r].island1);
        for(size_t r=0;r<nr;r++)if(relations[r].eligible){size_t root=(size_t)uf_find(&components,relations[r].island0);supported[root]=1;bad[root]|=relations[r].final_residual!=0;}
        for(size_t c=0;c<nc;c++)if(chart[c].valid && chart[c].phase_conflict)bad[uf_find(&components,(int32_t)piece_island[chart[c].piece])]=1;
    }
    double *minimum=ARENA_ALLOC(arena,(2*np+1)*sizeof *minimum);for(size_t g=0;g<=2*np;g++)minimum[g]=DBL_MAX;
    for(size_t c=0;c<nc;c++)if(chart[c].valid && !chart[c].phase_conflict){
        AroChart *q=chart+c;size_t island=piece_island[q->piece],group=(size_t)uf_find(&components,(int32_t)island);
        q->island=out->radial_solved?island:SIZE_MAX;
        q->relative=out->radial_solved && supported[group] && !bad[group];q->group=q->relative?group:np+q->piece;
        if(q->relative){q->qlo+=relative[island].lift;q->qhi+=relative[island].lift;}
        minimum[q->group]=fmin(minimum[q->group],q->qlo);
    }
    for(size_t c=0;c<nc;c++)if(AsmChart_registered(run->charts+c)){
        AroChart *q=chart+c;double area=run->charts[c].area3d;
        if(!q->valid || q->phase_conflict){out->unresolved_charts++;out->unresolved_area+=area;continue;}
        double origin=floor(minimum[q->group]);q->qlo-=origin;q->qhi-=origin;
        if(q->relative){out->relative_charts++;out->relative_area+=area;}else{out->within_charts++;out->within_area+=area;}
    }
    return 0;
}

static int aro_model(Arena_T arena,const AsmRun *run,const AsmAxis *axis,AroModel *out)
{return aro_model_mode(arena,run,axis,out,1);}

int AsmReadingOrder_ribbon_keys(const AsmRun *run,const AsmAxis *axis,AsmRibbonKey *keys)
{
    if(!run || !axis || !keys || !isfinite(run->layer_d_global) || !(run->layer_d_global>0))return -1;
    Arena_T arena=Arena_new();AroModel m={0};int rc=-1;
    if(aro_model_mode(arena,run,axis,&m,0))goto done;
    double *sum=ARENA_CALLOC(arena,2*(m.pieces+1),sizeof *sum);
    for(size_t c=0;c<run->n_charts;c++){
        const AroChart *q=m.chart+c;if(!q->valid || q->phase_conflict)continue;
        double w=run->charts[c].area3d;
        sum[2*q->piece]+=w*(q->radius/m.pitch-.5*(q->qlo+q->qhi));sum[2*q->piece+1]+=w;
    }
    for(size_t c=0;c<run->n_charts;c++){
        AsmRibbonKey *key=keys+c;memset(key,0,sizeof *key);
        if(!AsmChart_registered(run->charts+c))continue;
        const AroChart *q=m.chart+c;
        if(q->valid && !q->phase_conflict && sum[2*q->piece+1]>0){
            double offset=sum[2*q->piece]/sum[2*q->piece+1];
            key->turn_lo=q->qlo+offset;key->turn_hi=q->qhi+offset;key->source_supported=1;
            for(int d=0;d<2;d++)key->angular_gradient[d]=m.sense*q->angular_gradient[d];
        }else{
            double frame[3];AsmAxis_frame(axis,run->charts[c].centroid,frame);
            key->turn_lo=key->turn_hi=hypot(frame[1],frame[2])/m.pitch;
        }
        if(!isfinite(key->turn_lo) || !isfinite(key->turn_hi) ||
           !isfinite(key->angular_gradient[0]) || !isfinite(key->angular_gradient[1]))goto done;
    }
    rc=0;
done:
    Arena_dispose(&arena);return rc;
}

int AsmReadingOrder_write(const AsmRun *run,const AsmAxis *axis,const char *dir)
{
    if(!run || !axis || !dir)return -1;
    double started=ves_clock_sec();Arena_T arena=Arena_new();AroModel m={0};int rc=-1;FILE *fp=NULL;char path[2048];
    if(aro_model(arena,run,axis,&m))goto done;
    ves_mkdir(dir);snprintf(path,sizeof path,"%s/reading_order.json",dir);fp=fopen(path,"wb");if(!fp)goto done;
    fprintf(fp,"{\"schema\":\"native-physical-reading-order-v1\",\"source_modified\":false,\"identity_certified\":false,"
        "\"scope\":\"Relative turns only within each order_group. Source-cycle or radial contradictions prevent shared ordering. Atlas coordinates and missing material are unchanged.\","
        "\"pitch\":%.17g,\"sense\":%d,\"sense_evidence_pieces\":%zu,\"samples\":%zu,\"source_pieces\":%zu,\"source_edges\":%zu,\"source_cycle_conflicts\":%zu,"
        "\"radial_solved\":%s,\"radial_relations\":%zu,\"radial_satisfied\":%zu,\"radial_eligible\":%zu,\"radial_cycle_conflicts\":%zu,"
        "\"relative_charts\":%zu,\"within_piece_charts\":%zu,\"unresolved_charts\":%zu,\"local_pitch_charts\":%zu,"
        "\"relative_area\":%.17g,\"within_piece_area\":%.17g,\"unresolved_area\":%.17g,\"charts\":[\n",
        m.pitch,m.sense,m.sense_pieces,m.samples,m.pieces,m.source_edges,m.phase_conflicts,m.radial_solved?"true":"false",m.sync.relations,
        m.sync.satisfied_relations,m.sync.eligible_relations,m.radial_cycle_conflicts,m.relative_charts,m.within_charts,m.unresolved_charts,m.local_pitch_charts,
        m.relative_area,m.within_area,m.unresolved_area);
    size_t count=0;
    for(size_t c=0;c<run->n_charts;c++)if(AsmChart_registered(run->charts+c)){
        const AroChart *q=m.chart+c;int valid=q->valid && !q->phase_conflict;
        fprintf(fp,"%s{\"id\":%zu,\"status\":\"%s\",\"samples\":%zu,\"source_piece\":",count++?",\n":"",c,
                !valid?"unresolved":q->relative?"relative_group":"within_piece",q->count);
        if(q->valid)fprintf(fp,"%zu",q->piece);else fputs("null",fp);
        fputs(",\"radial_island\":",fp);if(q->island!=SIZE_MAX)fprintf(fp,"%zu",q->island);else fputs("null",fp);
        fputs(",\"order_group\":",fp);if(valid)fprintf(fp,"%zu",q->group);else fputs("null",fp);
        fputs(",\"turn\":",fp);if(valid)fprintf(fp,"[%.9g,%.9g]",q->qlo,q->qhi);else fputs("null",fp);
        fputs(",\"axial\":",fp);if(q->valid)fprintf(fp,"[%.9g,%.9g]",q->axial_lo,q->axial_hi);else fputs("null",fp);
        fprintf(fp,",\"radius\":%.9g,\"local_pitch\":%.9g,\"source_cycle_conflict\":%s}",q->radius,q->pitch,q->phase_conflict?"true":"false");
    }
    fprintf(fp,"\n],\"seconds\":%.9g}\n",ves_clock_sec()-started);rc=ferror(fp)?-1:0;
    if(fclose(fp))rc=-1;fp=NULL;if(rc)goto done;
    /* The local file viewer can load this wrapper without a web server. */
    FILE *json=fopen(path,"rb");if(!json){rc=-1;goto done;}
    snprintf(path,sizeof path,"%s/reading_order.js",dir);fp=fopen(path,"wb");
    if(!fp){fclose(json);rc=-1;goto done;}
    fputs("window.SHEET_READING_ORDER=",fp);char buffer[8192];size_t bytes=0;
    while((bytes=fread(buffer,1,sizeof buffer,json))!=0)if(fwrite(buffer,1,bytes,fp)!=bytes){rc=-1;break;}
    if(ferror(json))rc=-1;if(fclose(json))rc=-1;fputs(";\n",fp);
    if(ferror(fp))rc=-1;if(fclose(fp))rc=-1;fp=NULL;if(rc)goto done;
    snprintf(path,sizeof path,"%s/reading_relations.csv",dir);fp=fopen(path,"wb");if(!fp){rc=-1;goto done;}
    fputs("island_a,island_b,observations,cross_section_observations,seam_observations,target,mode_agreement,residual_median,eligible,relative_difference,cycle_residual\n",fp);
    for(size_t r=0;r<m.nrelations;r++){
        const AtlasWindingSyncRelation *e=m.relations+r;
        fprintf(fp,"%d,%d,%zu,%zu,%zu,%d,%.17g,%.17g,%d,%d,%d\n",e->island0,e->island1,e->observations,
                e->cross_section_observations,e->seam_observations,e->target_turn_correction,e->mode_agreement,e->residual_median,
                e->eligible,e->solved_turn_correction,e->final_residual);
    }
    if(ferror(fp))rc=-1;
    fprintf(stderr,"[reading order] %zu relative, %zu within-piece, %zu unresolved charts; %zu/%zu radial relations; %.3f s\n",
            m.relative_charts,m.within_charts,m.unresolved_charts,m.sync.satisfied_relations,m.sync.eligible_relations,ves_clock_sec()-started);
done:
    if(fp && fclose(fp))rc=-1;Arena_dispose(&arena);return rc;
}

#include "asm_reading_order_test.inc"
